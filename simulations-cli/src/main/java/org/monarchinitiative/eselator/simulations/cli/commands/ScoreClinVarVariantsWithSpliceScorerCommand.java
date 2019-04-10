package org.monarchinitiative.eselator.simulations.cli.commands;

import de.charite.compbio.jannovar.data.JannovarData;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.monarchinitiative.eselator.simulations.cli.CustomVariantEvaluation;
import org.monarchinitiative.exomiser.core.genome.VariantAnnotator;
import org.monarchinitiative.exomiser.core.genome.dao.splicing.*;
import org.monarchinitiative.exomiser.core.model.VariantAnnotation;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.phenopackets.schema.v1.Phenopacket;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static htsjdk.variant.vcf.VCFEncoder.formatVCFDouble;

@Component
public class ScoreClinVarVariantsWithSpliceScorerCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(ScoreClinVarVariantsWithSpliceScorerCommand.class);

    private static final Phenopacket DEFAULT_PHENOPACKET = Phenopacket.getDefaultInstance();

    private static final Pattern CLNSIG_REGEXP = Pattern.compile(".*Benign.*");

    private final SplicingParameters splicingParameters;

    private final SplicingInformationContentAnnotator splicingInformationContentAnnotator;

    private final JannovarData jannovarData;

    private final VariantAnnotator variantAnnotator;

    private final SplicingDao splicingDao;

    private final Path hexamerFilePath;

    private Path clinVarVcfPath, outputPath;

    public ScoreClinVarVariantsWithSpliceScorerCommand(JannovarData jannovarData,
                                                       VariantAnnotator variantAnnotator,
                                                       SplicingDao splicingDao, Path hexamerFilePath) {
        this.splicingParameters = splicingDao.getSplicingParameters();
        this.splicingInformationContentAnnotator = splicingDao.getIcAnnotator();
        this.jannovarData = jannovarData;
        this.variantAnnotator = variantAnnotator;
        this.splicingDao = splicingDao;
        this.hexamerFilePath = hexamerFilePath;
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

        // define predicate for selection of ClinVar variants for scoring
        Predicate<VariantContext> benignVariants = vc -> {
            String clnsig = vc.getAttributeAsString("CLNSIG", "Pathogenic");
            return CLNSIG_REGEXP.matcher(clnsig).matches();
        };

        // ----------------- SCORE VARIANTS & WRITE TO FILE -------------------
        LOGGER.info("Scoring variants");
        Collection<SpliceScorer> scorers = Scorers.getAllSpliceScorers(splicingParameters, splicingInformationContentAnnotator, hexamerFilePath);

        List<String> header = new ArrayList<>(Collections.singletonList("VARIANT"));
        List<String> scorerNames = scorers.stream().map(SpliceScorer::getName).sorted().collect(Collectors.toList());
        header.addAll(scorerNames);


        try (VCFFileReader reader = new VCFFileReader(clinVarVcfPath.toFile(), false);
             CSVPrinter printer = new CSVPrinter(Files.newBufferedWriter(outputPath),
                     CSVFormat.TDF.withHeader(header.toArray(new String[0])))
        ) {
            VCFHeader clinVarHeader = reader.getFileHeader();
            reader.iterator().stream()
                    .filter(benignVariants)
                    .forEach(vc -> {
                        if (vc.getNAlleles() < 2) {
                            // ALT allele is missing ('.'), nothing to be done here
                            return;
                        }
                        // ---------------------- SCORE ----------------------------------------------------------------
                        int chr = jannovarData.getRefDict().getContigNameToID().get(vc.getContig());

                        // clinvar vcf does not contain multiallelic variant records
                        VariantAnnotation annotation = variantAnnotator.annotate(vc.getContig(), vc.getStart(),
                                vc.getReference().getBaseString(), vc.getAlternateAllele(0).getBaseString());


                        VariantEvaluation ve = VariantEvaluation.builder(chr, vc.getStart(), vc.getReference().getBaseString(),
                                vc.getAlternateAllele(0).getBaseString())
                                .annotations(annotation.getTranscriptAnnotations())
                                .build();

                        final CustomVariantEvaluation.Builder builder = CustomVariantEvaluation.builder()
                                .setPhenopacket(DEFAULT_PHENOPACKET)
                                .setVariantEvaluation(ve)
                                .setVcfAlleleInfoField(MethodsFromHTSJDK.mapToVcfInfoField(vc.getAttributes(), clinVarHeader));

                        Optional<SplicingContext> sco = splicingDao.buildSplicingContext(ve);

                        if (!sco.isPresent()) {
                            // intergenic variants
                            return;
                        }

                        for (SpliceScorer scorer : scorers) {
                            try {
                                builder.putScore(scorer.getName(), scorer.getScoringFunction().apply(sco.get()));
                            } catch (IndexOutOfBoundsException e) {
                                LOGGER.warn("Error '{}'", sco.get(), e);
                                throw new RuntimeException(e);
                            }
                        }

                        CustomVariantEvaluation cve = builder.build();

                        // ---------------------- WRITE ----------------------------------------------------------------
                        List<Object> fields = new ArrayList<>();
                        // VARIANT
                        String variantRecord = String.format("%s:%d%s>%s", ve.getChromosomeName(), ve.getPosition(), ve.getRef(), ve.getAlt());
                        fields.add(variantRecord);

                        // SCORES
                        for (String name : scorerNames) {
                            fields.add(cve.getScoreMap().get(name));
                        }

                        // write the variant
                        try {
                            printer.printRecord(fields.toArray());
                        } catch (IOException e) {
                            LOGGER.warn("Unable to write variant '{}'", variantRecord, e);
                        }
                    });
        }
        LOGGER.info("くまくま━━━━━━ヽ（ ・(ｪ)・ ）ノ━━━━━━ !!!");
        LOGGER.info("                 Done!               ");
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
