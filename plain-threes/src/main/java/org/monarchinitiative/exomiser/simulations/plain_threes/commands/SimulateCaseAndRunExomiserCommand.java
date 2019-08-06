package org.monarchinitiative.exomiser.simulations.plain_threes.commands;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import de.charite.compbio.jannovar.mendel.SubModeOfInheritance;
import org.monarchinitiative.exomiser.core.Exomiser;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisMode;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.util.InheritanceModeOptions;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.model.GeneScore;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
import org.monarchinitiative.exomiser.core.writers.AnalysisResultsWriter;
import org.monarchinitiative.exomiser.core.writers.OutputFormat;
import org.monarchinitiative.exomiser.core.writers.OutputSettings;
import org.monarchinitiative.exomiser.simulations.plain_threes.Utils;
import org.monarchinitiative.exomiser.simulations.plain_threes.simulators.SingleVcfSimulator;
import org.monarchinitiative.exomiser.simulations.plain_threes.simulators.VcfSimulator;
import org.phenopackets.schema.v1.Phenopacket;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This runner implements command `--simulate-case-and-run-exomiser`.
 * <p>
 * Take a directory full of Phenopackets, simulate exome VCF based on a single VCF file, run Exomiser with and without
 * SPLICING score and write the results into given directory.
 */
@Component
public class SimulateCaseAndRunExomiserCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(SimulateCaseAndRunExomiserCommand.class);

    /**
     * For splicing aware analysis - SPLICING, REVEL and MVP only.
     */
    private static final Set<PathogenicitySource> PS_W_SPLICING = EnumSet.of(PathogenicitySource.REVEL, PathogenicitySource.MVP, PathogenicitySource.SPLICING);

    /**
     * For splicing agnostic analysis - REVEL and MVP only.
     */
    private static final Set<PathogenicitySource> PS_NOT_SPLICING = EnumSet.of(PathogenicitySource.REVEL, PathogenicitySource.MVP);

    /**
     * For frequency & inheritance filtering
     */
    private static final InheritanceModeOptions INHERITANCE_MODE_OPTIONS;

    private static final EnumSet<FrequencySource> FREQUENCY_SOURCES = EnumSet.of(
            FrequencySource.ESP_AFRICAN_AMERICAN, FrequencySource.ESP_ALL, FrequencySource.ESP_EUROPEAN_AMERICAN,
            FrequencySource.THOUSAND_GENOMES,
            FrequencySource.EXAC_AFRICAN_INC_AFRICAN_AMERICAN, FrequencySource.EXAC_AMERICAN, FrequencySource.EXAC_EAST_ASIAN, FrequencySource.EXAC_FINNISH,
            FrequencySource.EXAC_NON_FINNISH_EUROPEAN, FrequencySource.EXAC_SOUTH_ASIAN, FrequencySource.EXAC_OTHER,
            FrequencySource.UK10K, FrequencySource.TOPMED,
            FrequencySource.GNOMAD_E_AFR, FrequencySource.GNOMAD_E_AMR,
            // including the GNOMAD_E_ASJ and GNOMAD_G_ASJ reduces performance on 100K genomes so we disable it by default
            //FrequencySource.GNOMAD_E_ASJ,
            FrequencySource.GNOMAD_E_EAS, FrequencySource.GNOMAD_E_FIN,
            FrequencySource.GNOMAD_E_NFE, FrequencySource.GNOMAD_E_OTH, FrequencySource.GNOMAD_E_SAS,
            FrequencySource.GNOMAD_G_AFR, FrequencySource.GNOMAD_G_AMR,
            //FrequencySource.GNOMAD_G_ASJ,
            FrequencySource.GNOMAD_G_EAS, FrequencySource.GNOMAD_G_FIN,
            FrequencySource.GNOMAD_G_NFE, FrequencySource.GNOMAD_G_OTH, FrequencySource.GNOMAD_G_SAS);

    private static final EnumSet<VariantEffect> NON_CODING_EFFECTS = EnumSet.of(
//            VariantEffect.FIVE_PRIME_UTR_EXON_VARIANT,
//            VariantEffect.FIVE_PRIME_UTR_INTRON_VARIANT,
//            VariantEffect.THREE_PRIME_UTR_EXON_VARIANT,
//            VariantEffect.THREE_PRIME_UTR_INTRON_VARIANT,
            VariantEffect.NON_CODING_TRANSCRIPT_EXON_VARIANT,
            VariantEffect.NON_CODING_TRANSCRIPT_INTRON_VARIANT,
//            VariantEffect.CODING_TRANSCRIPT_INTRON_VARIANT,
            VariantEffect.UPSTREAM_GENE_VARIANT,
            VariantEffect.DOWNSTREAM_GENE_VARIANT,
            VariantEffect.INTERGENIC_VARIANT,
            VariantEffect.REGULATORY_REGION_VARIANT
    );

    private static final Set<OutputFormat> OUTPUT_FORMATS = EnumSet.of(OutputFormat.HTML, OutputFormat.TSV_VARIANT, OutputFormat.VCF);

    static {
        Map<SubModeOfInheritance, Float> inheritanceModeFrequencyCutoffs = new EnumMap<>(SubModeOfInheritance.class);
        // all frequencies are in percentage values
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.AUTOSOMAL_DOMINANT, 0.1f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.AUTOSOMAL_RECESSIVE_COMP_HET, 2.0f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.AUTOSOMAL_RECESSIVE_HOM_ALT, 0.1f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.X_DOMINANT, 0.1f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.X_RECESSIVE_COMP_HET, 2.0f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.X_RECESSIVE_HOM_ALT, 0.1f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.MITOCHONDRIAL, 0.2f);
        INHERITANCE_MODE_OPTIONS = InheritanceModeOptions.of(inheritanceModeFrequencyCutoffs);
    }

    // ------------------------------      DEPENDENCIES      ------------------------------------------
    private final Exomiser exomiser;


    // ------------------------------        CLI ARGS        ------------------------------------------

    /**
     * List of paths pointing to phenopackets that should be evaluated.
     */
    private final List<Path> phenopacketPaths = new ArrayList<>();

    /**
     * Path to VCF with background variants used to simulate an exome.
     */
    private Path templateVcfPath;

    /**
     * Path to directory where output will be directed.
     */
    private Path outputPath;


    public SimulateCaseAndRunExomiserCommand(Exomiser exomiser) {
        this.exomiser = exomiser;
    }

    private static SimpleResults evaluateResults(Phenopacket pp, AnalysisResults splicingAgnosticResults, AnalysisResults splicingAwareResults) {
        String geneSymbol = pp.getGenes(0).getSymbol();

        int agnosticRank = -1;
        final List<GeneScore> agnosticScores = splicingAgnosticResults.getGeneScores();
        for (int i = 0; i < agnosticScores.size(); i++) {
            final GeneScore gs = agnosticScores.get(i);
            if (gs.getGeneIdentifier().getHgncSymbol().equals(geneSymbol)) {
                agnosticRank = i + 1; // if i=0, then the gene was in fact the gene #1
                break;
            }
        }

        int awareRank = -1;
        final List<GeneScore> awareScores = splicingAwareResults.getGeneScores();
        for (int i = 0; i < awareScores.size(); i++) {
            final GeneScore gs = awareScores.get(i);
            if (gs.getGeneIdentifier().getHgncSymbol().equals(geneSymbol)) {
                awareRank = i + 1;
                break;
            }
        }

        return SimpleResults.builder()
                .setCaseName(pp.getId())
                .setGeneSymbol(geneSymbol)
                .setSplicingAgnosticRank(agnosticRank)
                .setSplicingAwareRank(awareRank)
                .setSplicingPathomechanism(Utils.getSplicingPathomechanism(pp.getVariantsList()))
                .build();
    }

    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("simulate-case-and-run-exomiser")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }

        VcfSimulator simulator = new SingleVcfSimulator(templateVcfPath);

        Path ranksPath = outputPath.resolve("ranks.tsv");
        try (BufferedWriter resultWriter = Files.newBufferedWriter(ranksPath)) {
            // write header of the ranks file
            String delimiter = "\t";
            resultWriter.write(String.join(delimiter, Arrays.asList("CASE", "WITH_SPLICING", "WITHOUT_SPLICING", "PATHOMECHANISM")));
            resultWriter.newLine();

            // -----------------------    FOR EACH PHENOPACKET    --------------------------------------
            for (Path phenopacketPath : phenopacketPaths) {
                LOGGER.info("Reading phenopacket from '{}'", phenopacketPath);
                Phenopacket pp = Utils.readPhenopacket(phenopacketPath);
                if (pp.getSubject().getId().isEmpty()) {
                    LOGGER.error("Phenopacket subject's ID must not be empty. Unable to continue");
                    continue;
                }

                if (pp.getGenesCount() != 1) {
                    LOGGER.error("Phenopackets used for simulation MUST have exactly one gene");
                    continue;
                }

                if (pp.getDiseasesCount() != 1) {
                    LOGGER.error("Phenopackets used for simulation MUST have exactly one disease");
                    continue;
                }


                // -----------------------    CREATE THE SIMULATED VCF FILE    -------------------------
                LOGGER.info("Creating simulated VCF file");
                Path vcfPath = simulator.simulateVcfWithPhenopacket(pp);

                String ppFileName = phenopacketPath.toFile().getName();

                // Exomiser results for given phenopacket will be written here
                Path phenopacketOutputDir = Files.createDirectories(this.outputPath.resolve(ppFileName));

                // -------------------------------------------------------------------------------------
                //
                // First run the analysis without splicing scores
                //
                LOGGER.info("\n\n\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A" +
                        "   Creating splicing-agnostic analysis   " +
                        "\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\u266B\u266C\u266A\n");
                String sampleName = pp.getSubject().getId().replaceAll("\\s+", "_");

                List<String> phenotypesAsHpoStrings = Utils.getPresentPhenotypesAsHpoStrings(pp);

                Analysis splicingAgnosticAnalysis = exomiser.getAnalysisBuilder()
                        .genomeAssembly(GenomeAssembly.HG19)
                        .vcfPath(vcfPath)
                        .probandSampleName(sampleName)
                        .hpoIds(phenotypesAsHpoStrings)
                        .analysisMode(AnalysisMode.PASS_ONLY)
                        .inheritanceModes(INHERITANCE_MODE_OPTIONS)
                        .frequencySources(FREQUENCY_SOURCES)
                        .pathogenicitySources(PS_NOT_SPLICING) // all the pathogenicity sources except SPLICING & TEST
                        // adds an mask for removing non-coding variants
                        .addQualityFilter(200)
                        .addVariantEffectFilter(NON_CODING_EFFECTS)
                        .addFailedVariantFilter()
                        // frequency filter max will be automatically derived from the inheritance mode options
                        .addFrequencyFilter()
                        .addPathogenicityFilter(true)
                        .addInheritanceFilter()
                        .addOmimPrioritiser()
                        .addHiPhivePrioritiser()
                        .build();

                AnalysisResults splicingAgnosticResults = exomiser.run(splicingAgnosticAnalysis);

                OutputSettings agnosticSettings = OutputSettings.builder()
                        .outputFormats(OUTPUT_FORMATS)
                        .outputPrefix(phenopacketOutputDir.resolve(ppFileName + "_NO").toString())
                        .build();
                AnalysisResultsWriter.writeToFile(splicingAgnosticAnalysis, splicingAgnosticResults, agnosticSettings);

                //
                // Then run the analysis with splicing pathogenicity scores
                LOGGER.info("\n\n\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708\u2708" +
                        "   Creating splicing-aware analysis   " +
                        "\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\u2600\n");
                Analysis splicingAwareAnalysis = exomiser.getAnalysisBuilder()
                        .genomeAssembly(GenomeAssembly.HG19)
                        .vcfPath(vcfPath)
                        .probandSampleName(sampleName)
                        .hpoIds(phenotypesAsHpoStrings)
                        .analysisMode(AnalysisMode.PASS_ONLY)
                        .inheritanceModes(INHERITANCE_MODE_OPTIONS)
                        .frequencySources(FrequencySource.ALL_EXTERNAL_FREQ_SOURCES)
                        .pathogenicitySources(PS_W_SPLICING) // all the pathogenicity sources except TEST
                        // adds an mask for removing non-coding variants
                        .addQualityFilter(200)
                        .addVariantEffectFilter(NON_CODING_EFFECTS)
                        .addFailedVariantFilter()
                        // frequency filter max will be automatically derived from the inheritance mode options
                        .addFrequencyFilter()
                        .addPathogenicityFilter(true)
                        .addInheritanceFilter()
                        .addOmimPrioritiser()
                        .addHiPhivePrioritiser()
                        .build();

                AnalysisResults splicingAwareResults = exomiser.run(splicingAwareAnalysis);

                OutputSettings awareSettings = OutputSettings.builder()
                        .outputFormats(OUTPUT_FORMATS)
                        .outputPrefix(phenopacketOutputDir.resolve(ppFileName + "_YES").toString())
                        .build();
                AnalysisResultsWriter.writeToFile(splicingAwareAnalysis, splicingAwareResults, awareSettings);


                //
                // write ranks/evaluation of the analyses
                SimpleResults sr = evaluateResults(pp, splicingAgnosticResults, splicingAwareResults);
                String resultLine = sr.getCaseName() + delimiter
                        + sr.getRankWithSplicing() + delimiter
                        + sr.getRankNormal() + delimiter
                        + sr.getSplicingPathomechanism();
                resultWriter.write(resultLine);
                resultWriter.newLine();
                resultWriter.flush();
            }
        }

        LOGGER.info("くまくま━━━━━━ヽ（ ・(ｪ)・ ）ノ━━━━━━ !!!");
        LOGGER.info("                 Done!               ");
    }


    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopackets
        if (args.containsOption("pp-dir")) {
            String ppDirString = args.getOptionValues("pp-dir").get(0);
            Path ppDirPath = Paths.get(ppDirString);
            File[] jsonFiles = ppDirPath.toFile().listFiles(f -> f.getName().endsWith(".json"));
            if (jsonFiles != null) {
                Arrays.stream(jsonFiles).map(File::toPath).forEach(phenopacketPaths::add);
            }
        }

        List<String> pps = args.getOptionValues("pp");
        if (pps != null) {
            phenopacketPaths.addAll(pps.stream().map(Paths::get).collect(Collectors.toList()));
        }


        //
        if (!args.containsOption("vcf")) {
            LOGGER.error("Missing 'vcf' argument");
            return false;
        }
        templateVcfPath = Paths.get(args.getOptionValues("vcf").get(0));

        // Output directory path - where to write all the results
        if (!args.containsOption("output-exomiser")) {
            LOGGER.error("Missing '--output-exomiser' argument");
            return false;
        }
        outputPath = Paths.get(args.getOptionValues("output-exomiser").get(0));

        return true;
    }

    private static class SimpleResults {

        /**
         * Name of the case represented by phenopacket.
         */
        private final String caseName;

        /**
         * String with HGNC symbol of the causal gene, e.g. 'GCK1'.
         */
        private final String geneSymbol;

        /**
         * Rank of the causal gene created by either with splicing aware exomiser analysis or splicing agnostic analysis.
         */
        private final int rankWithSplicing, rankNormal;

        /**
         * Background information about the splicing variant, e.g. `splicing|3ss|disrupted`
         */
        private final String splicingPathomechanism;

        private SimpleResults(Builder builder) {
            caseName = builder.caseName;
            geneSymbol = builder.geneSymbol;
            rankWithSplicing = builder.splicingAwareRank;
            rankNormal = builder.splicingAgnosticRank;
            splicingPathomechanism = builder.splicingPathomechanism;
        }

        private static SimpleResults.Builder builder() {
            return new SimpleResults.Builder();
        }

        @Override
        public String toString() {
            return "SimpleResults{" +
                    "caseName='" + caseName + '\'' +
                    ", geneSymbol='" + geneSymbol + '\'' +
                    ", rankWithSplicing=" + rankWithSplicing +
                    ", rankNormal=" + rankNormal +
                    ", splicingPathomechanism='" + splicingPathomechanism + '\'' +
                    '}';
        }

        public String getCaseName() {
            return caseName;
        }

        public String getGeneSymbol() {
            return geneSymbol;
        }

        public int getRankWithSplicing() {
            return rankWithSplicing;
        }

        public int getRankNormal() {
            return rankNormal;
        }

        public String getSplicingPathomechanism() {
            return splicingPathomechanism;
        }

        public static final class Builder {

            private String caseName;

            private String geneSymbol;

            private int splicingAwareRank;

            private int splicingAgnosticRank;

            private String splicingPathomechanism;

            private Builder() {
            }

            public Builder setCaseName(String caseName) {
                this.caseName = caseName;
                return this;
            }

            public Builder setGeneSymbol(String geneSymbol) {
                this.geneSymbol = geneSymbol;
                return this;
            }

            public Builder setSplicingAwareRank(int splicingAwareRank) {
                this.splicingAwareRank = splicingAwareRank;
                return this;
            }

            public Builder setSplicingAgnosticRank(int splicingAgnosticRank) {
                this.splicingAgnosticRank = splicingAgnosticRank;
                return this;
            }

            public Builder setSplicingPathomechanism(String splicingPathomechanism) {
                this.splicingPathomechanism = splicingPathomechanism;
                return this;
            }

            public SimpleResults build() {
                return new SimpleResults(this);
            }
        }
    }
}
