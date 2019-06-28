package org.monarchinitiative.exomiser.simulations.plain_threes.commands;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
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
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Predicate;
import java.util.regex.Pattern;

import static htsjdk.variant.vcf.VCFEncoder.formatVCFDouble;

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

        List<String> header = new ArrayList<>(Collections.singletonList("VARIANT"));
        scoringStrategies.forEach(ss -> header.add(ss.toString()));


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
                List<SplicingTranscript> transcripts = splicingTranscriptSource.fetchTranscripts(varCoordinates.getContig(), varCoordinates.getBegin(), varCoordinates.getEnd());

                if (transcripts.isEmpty()) {
                    LOGGER.warn("No transcript overlaps with variant {}", splv);
                    continue;
                }

                SplicingTranscript longestOp = transcripts.stream()
                        .max(Comparator.comparing(SplicingTranscript::getTxLength))
                        .get();


                // fetch nucleotide sequence neighboring the variant
                SequenceInterval sequenceInterval = genomeSequenceAccessor.fetchSequence(longestOp.getContig(),
                        longestOp.getTxBegin() - 50, longestOp.getTxEnd() + 50,
                        longestOp.getStrand());


                // evaluate variant against all the transcripts & write out the scores
                for (SplicingTranscript transcript : transcripts) {
                    // evaluate
                    SplicingPathogenicityData evaluation = splicingEvaluator.evaluate(splv, transcript, sequenceInterval);

                    // write out
                    StringBuilder builder = new StringBuilder()
                            .append(String.format("%s:%d %s>%s", splv.getContig(), splv.getPos(), splv.getRef(), splv.getAlt())).append(DELIMITER)
                            .append(transcript.getAccessionId()).append(DELIMITER)
                            .append(evaluation.getMaxScore()); // no delimiter here!

                    // write out all the scores
                    for (ScoringStrategy strategy : scoringStrategies) {
                        builder.append(DELIMITER).append(evaluation.getScoresMap().getOrDefault(strategy, Double.NaN));
                    }

                    writer.write(builder.toString());
                    writer.newLine();
                }
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
        if (!args.containsOption("output")) {
            LOGGER.warn("Missing '--output' argument");
            return false;
        }
        outputPath = Paths.get(args.getOptionValues("output").get(0));

        // Strict flag
        strict = args.containsOption("strict");

        return true;
    }


    /**
     * Methods related to variant context copied from HtsJDK {@link VCFEncoder} class.
     */
    private static class MethodsFromHTSJDK {


        private static String mapToVcfInfoField(Map<String, Object> attributes, VCFHeader header) {
            final Map<String, String> infoFields = new TreeMap<>();
            for (final Map.Entry<String, Object> field : attributes.entrySet()) {
//            We expect the ClinVar VCF to be well-formatted
//            if (!header.hasInfoLine(field.getKey()))
//                fieldIsMissingFromHeaderError(context, field.getKey(), "INFO");

                final String outputValue = formatVCFField(field.getValue());
                if (outputValue != null) infoFields.put(field.getKey(), outputValue);
            }
            return createInfoString(infoFields, header);
        }

        private static String formatVCFField(final Object val) {
            final String result;
            if (val == null)
                result = VCFConstants.MISSING_VALUE_v4;
            else if (val instanceof Double)
                result = formatVCFDouble((Double) val);
            else if (val instanceof Boolean)
                result = (Boolean) val ? "" : null; // empty string for true, null for false
            else if (val instanceof List) {
                result = formatVCFField(((List) val).toArray());
            } else if (val.getClass().isArray()) {
                final int length = Array.getLength(val);
                if (length == 0)
                    return formatVCFField(null);
                final StringBuilder sb = new StringBuilder(formatVCFField(Array.get(val, 0)));
                for (int i = 1; i < length; i++) {
                    sb.append(',');
                    sb.append(formatVCFField(Array.get(val, i)));
                }
                result = sb.toString();
            } else
                result = val.toString();

            return result;
        }


        /*
         * Create the info string; assumes that no values are null
         */
        private static String createInfoString(final Map<String, String> infoFields, final VCFHeader header) {
            if (infoFields.isEmpty()) {
                return VCFConstants.EMPTY_INFO_FIELD;
            }

            StringBuilder builder = new StringBuilder();

            boolean isFirst = true;
            for (final Map.Entry<String, String> entry : infoFields.entrySet()) {
                if (isFirst) isFirst = false;
                else builder.append(VCFConstants.INFO_FIELD_SEPARATOR);

                builder.append(entry.getKey());

                if (!entry.getValue().equals("")) {
                    final VCFInfoHeaderLine metaData = header.getInfoHeaderLine(entry.getKey());
                    if (metaData == null || metaData.getCountType() != VCFHeaderLineCount.INTEGER || metaData.getCount() != 0) {
                        builder.append('=');
                        builder.append(entry.getValue());
                    }
                }
            }
            return builder.toString();
        }
    }
}
