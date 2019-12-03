package org.monarchinitiative.threes.benchmarks.cmd;

import de.charite.compbio.jannovar.data.ReferenceDictionary;
import de.charite.compbio.jannovar.reference.*;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import net.sourceforge.argparse4j.inf.Namespace;
import net.sourceforge.argparse4j.inf.Subparser;
import net.sourceforge.argparse4j.inf.Subparsers;
import org.monarchinitiative.threes.benchmarks.Utils;
import org.monarchinitiative.threes.core.data.SplicingTranscriptSource;
import org.monarchinitiative.threes.core.model.SplicingTranscript;
import org.monarchinitiative.threes.core.scoring.SplicingAnnotator;
import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;
import org.monarchinitiative.threes.core.scoring.sparse.ScoringStrategy;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;
import xyz.ielis.hyperutil.reference.fasta.GenomeSequenceAccessor;
import xyz.ielis.hyperutil.reference.fasta.SequenceInterval;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
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
public class ClinvarScorerCommand extends Command {

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

    private final SplicingAnnotator splicingAnnotator;

    public ClinvarScorerCommand(GenomeSequenceAccessor genomeSequenceAccessor,
                                SplicingTranscriptSource splicingTranscriptSource,
                                SplicingAnnotator splicingAnnotator) {
        this.genomeSequenceAccessor = genomeSequenceAccessor;
        this.splicingTranscriptSource = splicingTranscriptSource;
        this.splicingAnnotator = splicingAnnotator;
    }


    public static void setupSubparsers(Subparsers subparsers) {
        Subparser ic = subparsers.addParser("clinvar-scorer")
                .setDefault("cmd", "clinvar-scorer")
                .help("run benchmarks with ClinVar VCF");

        ic.addArgument("clinvar-vcf")
                .help("path to ClinVar VCF file");

        ic.addArgument("output-clinvar")
                .help("where to write the output file");

        ic.addArgument("-s", "--strict")
                .type(Boolean.class)
                .setDefault(false)
                .help("run in strict mode");
    }


    @Override
    public void run(Namespace namespace) throws CommandException {
        // ----------------------       CLI ARGS       ------------------------------------------------------------------
        Path clinVarVcfPath = Paths.get(namespace.getString("clinvar_vcf"));
        Path outputPath = Paths.get(namespace.getString("output_clinvar"));
        boolean strict = namespace.getBoolean("strict");

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

        final ReferenceDictionary rd = genomeSequenceAccessor.getReferenceDictionary();
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
                if (!rd.getContigNameToID().containsKey(vc.getContig())) {
                    LOGGER.warn("Unknown chromosome `{}` for variant {}:{}{}>{}", vc.getContig(), vc.getContig(), vc.getStart(), vc.getReference().getBaseString(), vc.getAlternateAlleles().get(0).getBaseString());
                    continue;
                }
                for (Allele altAllele : vc.getAlternateAlleles()) {
                    final GenomeVariant gv = new GenomeVariant(new GenomePosition(rd, Strand.FWD, rd.getContigNameToID().get(vc.getContig()), vc.getStart(), PositionType.ONE_BASED),
                            vc.getReference().getBaseString(), altAllele.getBaseString());
                    LOGGER.info("Evaluating `{}`", gv);


                    // fetch all the transcripts overlapping with variant's position
                    // retain the curated transcripts (accession id starts with 'NM_')
                    List<SplicingTranscript> curatedRefseqTranscripts = splicingTranscriptSource.fetchTranscripts(gv.getChrName(), gv.getGenomeInterval().getBeginPos(), gv.getGenomeInterval().getEndPos(), rd).stream()
                            .filter(tx -> tx.getAccessionId().startsWith("NM_"))
                            .collect(Collectors.toList());

                    if (curatedRefseqTranscripts.isEmpty()) {
                        LOGGER.warn("No curated transcript overlaps with variant `{}`", gv);
                        continue;
                    }


                    SplicingTranscript transcript = curatedRefseqTranscripts.stream()
                            .min(Utils.transcriptPriorityComparator())
                            .get();


                    // fetch nucleotide sequence neighboring the variant
                    final GenomeInterval txQueryRegion = transcript.getTxRegionCoordinates().withMorePadding(50, 50);
                    Optional<SequenceInterval> sequenceIntervalOpt = genomeSequenceAccessor.fetchSequence(txQueryRegion);
                    if (!sequenceIntervalOpt.isPresent()) {
                        LOGGER.warn("Unable to fetch sequence `{}` for transcript `{}-{}`", txQueryRegion, transcript.getAccessionId(), transcript.getTxRegionCoordinates());
                        continue;
                    }


                    // --- EVALUATE VARIANT AGAINST THE TRANSCRIPT & WRITE OUT THE SCORES ---
                    // evaluate
                    SplicingPathogenicityData evaluation = splicingAnnotator.evaluate(gv, transcript, sequenceIntervalOpt.get());

                    // write out
                    StringBuilder builder = new StringBuilder()
                            .append(String.format("%s:%d %s>%s", gv.getChrName(), gv.getPos() + 1, gv.getRef(), gv.getAlt())).append(DELIMITER)
                            .append(transcript.getAccessionId()).append(DELIMITER)
                            .append(evaluation.getMaxScore()); // no delimiter here!

                    // write out all the scores
                    for (ScoringStrategy ss : scoringStrategies) {
                        builder.append(DELIMITER).append(evaluation.getScoresMap().getOrDefault(ss, Double.NaN));
                    }

                    writer.write(builder.toString());
                    writer.newLine();
                }
            }

            LOGGER.info("くまくま━━━━━━ヽ（ ・(ｪ)・ ）ノ━━━━━━ !!!");
            LOGGER.info("                 Done!               ");

        } catch (IOException e) {
            LOGGER.error("Error: ", e);
            throw new CommandException(e);
        }
    }
}
