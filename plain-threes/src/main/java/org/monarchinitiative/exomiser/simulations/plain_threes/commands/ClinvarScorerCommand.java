package org.monarchinitiative.exomiser.simulations.plain_threes.commands;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.monarchinitiative.exomiser.simulations.plain_threes.Utils;
import org.monarchinitiative.threes.core.data.SplicingTranscriptSource;
import org.monarchinitiative.threes.core.model.GenomeCoordinates;
import org.monarchinitiative.threes.core.model.SequenceInterval;
import org.monarchinitiative.threes.core.model.SplicingTranscript;
import org.monarchinitiative.threes.core.model.SplicingVariant;
import org.monarchinitiative.threes.core.reference.fasta.GenomeSequenceAccessor;
import org.monarchinitiative.threes.core.scoring.ScoringStrategy;
import org.monarchinitiative.threes.core.scoring.SplicingEvaluator;
import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;


/**
 * This command runs the `--clinvar-scorer` command.
 * <p>
 * The command takes ClinVar VCF file, selects variants with benign or likely benign clinical significance (see `--strict`
 * flag) and scores variants using all splicing strategies.<br>
 * The results are written into a tsv file.
 * </p>
 */
@Component
public class ClinvarScorerCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(ClinvarScorerCommand.class);

    private static final String DELIMITER = "\t";

    /**
     * Matches 'Benign'
     */
    private static final Pattern BENIGN_CLNSIG = Pattern.compile(".*Benign.*");

    /**
     * Matches 'Likely_benign'
     */
    private static final Pattern LIKELY_BENIGN_CLNSIG = Pattern.compile(".*Likely_benign.*");

    // ----------------------      DEPENDENCIES    ------------------------------------------------------------------
    private final GenomeSequenceAccessor genomeSequenceAccessor;

    private final SplicingTranscriptSource splicingTranscriptSource;

    private final SplicingEvaluator splicingEvaluator;

    // ----------------------       CLI ARGS       ------------------------------------------------------------------
    private Path clinVarVcfPath;

    private Path outputPath;

    private boolean strict;

    public ClinvarScorerCommand(GenomeSequenceAccessor genomeSequenceAccessor, SplicingTranscriptSource splicingTranscriptSource, SplicingEvaluator splicingEvaluator) {
        this.genomeSequenceAccessor = genomeSequenceAccessor;
        this.splicingTranscriptSource = splicingTranscriptSource;
        this.splicingEvaluator = splicingEvaluator;
    }


    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("clinvar-scorer")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }

        // Analyze only Benign variants when `--strict` flag is present.
        // Analyze Benign & Likely benign variants without the `--strict` flag.
        Predicate<VariantContext> clnsigVariantMatcher = strict
                ? vc -> {
            String clnsig = vc.getAttributeAsString("CLNSIG", "Crap");
            return BENIGN_CLNSIG.matcher(clnsig).matches();
        }
                : vc -> {
            String clnsig = vc.getAttributeAsString("CLNSIG", "Crap");
            return BENIGN_CLNSIG.matcher(clnsig).matches() || LIKELY_BENIGN_CLNSIG.matcher(clnsig).matches();
        };


        List<ScoringStrategy> scoringStrategies = Arrays.asList(
                // donor
                ScoringStrategy.CANONICAL_DONOR, ScoringStrategy.CRYPTIC_DONOR, ScoringStrategy.CRYPTIC_DONOR_IN_CANONICAL_POSITION,
                // acceptor
                ScoringStrategy.CANONICAL_ACCEPTOR, ScoringStrategy.CRYPTIC_ACCEPTOR, ScoringStrategy.CRYPTIC_ACCEPTOR_IN_CANONICAL_POSITION,
                // ESE/ESS
                ScoringStrategy.SMS
        );

        // ----------------- SCORE VARIANTS & WRITE TO FILE -------------------
        LOGGER.info("Scoring variants");

        List<String> header = new ArrayList<>(Arrays.asList("VARIANT", "TX_ACC_ID", "MAX_SCORE"));
        for (ScoringStrategy ss : scoringStrategies) {
            header.add(ss.toString());
        }

        try (VCFFileReader reader = new VCFFileReader(clinVarVcfPath.toFile(), false);
             BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
            // write header
            writer.write(String.join(DELIMITER, header));
            writer.newLine();

            // for each variant in ClinVar
            for (VariantContext vc : reader) {
                // only process Benign variants
                if (!clnsigVariantMatcher.test(vc)) {
                    continue;
                }

                if (vc.getNAlleles() < 2) {
                    // ALT allele is missing ('.'), nothing to be done here
                    continue;
                }


                // make variant proper for splicing analysis
                GenomeCoordinates varCoordinates = GenomeCoordinates.newBuilder()
                        .setContig(vc.getContig())
                        .setBegin(vc.getStart() - 1)
                        .setEnd(vc.getStart() + vc.getReference().getBaseString().length() - 1)
                        .setStrand(true)
                        .build();

                SplicingVariant splv = SplicingVariant.newBuilder()
                        .setCoordinates(varCoordinates)
                        .setRef(vc.getReference().getBaseString())
                        .setAlt(vc.getAlternateAllele(0).getBaseString())
                        .build();
                LOGGER.info("Evaluating {}", splv);

                // fetch all the transcripts overlapping with variant's position
                // retain the curated transcripts (accession id starts with 'NM_')
                List<SplicingTranscript> curatedTranscripts = splicingTranscriptSource.fetchTranscripts(varCoordinates.getContig(), varCoordinates.getBegin(), varCoordinates.getEnd()).stream()
                        .filter(tx -> tx.getAccessionId().startsWith("NM_"))
                        .collect(Collectors.toList());

                if (curatedTranscripts.isEmpty()) {
                    LOGGER.warn("No curated transcript overlaps with variant {}", splv);
                    continue;
                }


                SplicingTranscript transcript = curatedTranscripts.stream()
                        .min(Utils.transcriptPriorityComparator())
                        .get();


                // fetch nucleotide sequence neighboring the variant
                SequenceInterval sequenceInterval = genomeSequenceAccessor.fetchSequence(transcript.getContig(),
                        transcript.getTxBegin() - 50, transcript.getTxEnd() + 50,
                        transcript.getStrand());


                // --- EVALUATE VARIANT AGAINST THE TRANSCRIPT & WRITE OUT THE SCORES ---
                // evaluate
                SplicingPathogenicityData evaluation = splicingEvaluator.evaluate(splv, transcript, sequenceInterval);

                // write out
                StringBuilder builder = new StringBuilder()
                        .append(String.format("%s:%d %s>%s", splv.getContig(), splv.getPos(), splv.getRef(), splv.getAlt())).append(DELIMITER)
                        .append(transcript.getAccessionId()).append(DELIMITER)
                        .append(evaluation.getMaxScore()); // no delimiter here!

                // write out all the scores
                for (ScoringStrategy ss : scoringStrategies) {
                    builder.append(DELIMITER).append(evaluation.getScoresMap().getOrDefault(ss, Double.NaN));
                }

                writer.write(builder.toString());
                writer.newLine();
            }

            LOGGER.info("くまくま━━━━━━ヽ（ ・(ｪ)・ ）ノ━━━━━━ !!!");
            LOGGER.info("                 Done!               ");

        }
    }


    private boolean parseCliArgs(ApplicationArguments args) {
        // Path to Clinvar vcf
        if (!args.containsOption("clinvar-vcf")) {
            LOGGER.warn("Missing '--clinvar-vcf' argument");
            return false;
        }
        clinVarVcfPath = Paths.get(args.getOptionValues("clinvar-vcf").get(0));


        // Output file path - results
        if (!args.containsOption("output-clinvar")) {
            LOGGER.warn("Missing '--output-clinvar' argument");
            return false;
        }
        outputPath = Paths.get(args.getOptionValues("output-clinvar").get(0));

        // Strict flag
        strict = args.containsOption("strict");

        return true;
    }

}
