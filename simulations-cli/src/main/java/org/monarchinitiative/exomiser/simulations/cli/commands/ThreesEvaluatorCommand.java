package org.monarchinitiative.exomiser.simulations.cli.commands;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.monarchinitiative.exomiser.core.Exomiser;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisMode;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.util.InheritanceModeOptions;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.model.GeneScore;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
import org.monarchinitiative.exomiser.core.prioritisers.model.Disease;
import org.monarchinitiative.exomiser.core.prioritisers.model.InheritanceMode;
import org.monarchinitiative.exomiser.core.prioritisers.service.PriorityService;
import org.monarchinitiative.exomiser.core.writers.AnalysisResultsWriter;
import org.monarchinitiative.exomiser.core.writers.OutputFormat;
import org.monarchinitiative.exomiser.core.writers.OutputSettings;
import org.monarchinitiative.exomiser.simulations.cli.Utils;
import org.monarchinitiative.exomiser.simulations.cli.simulators.SingleVcfSimulator;
import org.monarchinitiative.exomiser.simulations.cli.simulators.VcfSimulator;
import org.phenopackets.schema.v1.Phenopacket;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 *
 */
//@Component
@Deprecated // in favor of commands in the `plain-threes` modul

public class ThreesEvaluatorCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(ThreesEvaluatorCommand.class);

    private static final Set<PathogenicitySource> PS_W_SPLICING = Arrays.stream(PathogenicitySource.values())
            .filter(ps -> !ps.equals(PathogenicitySource.TEST))
            .collect(Collectors.toSet());

    private static final Set<PathogenicitySource> PS_NOT_SPLICING = PS_W_SPLICING.stream()
            .filter(ps -> !ps.equals(PathogenicitySource.SPLICING))
            .collect(Collectors.toSet());


    private static final float FREQ_CUTOFF = 1.0F;

    // ------------------------------      DEPENDENCIES      ------------------------------------------
    private final Exomiser exomiser;

    private final PriorityService priorityService;


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
     * Set to true if you want to evaluate with appropriate inheritance mode only. The inheritance mode is derived from
     * disease(s) associated with the causal gene.
     */
    private boolean inheritanceModeAware;


    /**
     * Path to directory where output will be directed.
     */
    private Path outputPath;


    public ThreesEvaluatorCommand(Exomiser exomiser, PriorityService priorityService) {
        this.exomiser = exomiser;
        this.priorityService = priorityService;
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

    private static void writeResults(List<SimpleResults> results, Path outputDirPath) throws IOException {
        String delimiter = "\t";
        try (final BufferedWriter writer = Files.newBufferedWriter(outputDirPath.resolve("ranks.tsv"))) {
            writer.write(String.join(delimiter, Arrays.asList("CASE", "WITH_SPLICING", "WITHOUT_SPLICING", "PATHOMECHANISM")));
            writer.newLine();
            for (SimpleResults result : results) {
                writer.write(result.getCaseName() + delimiter
                        + result.getRankWithSplicing() + delimiter
                        + result.getRankNormal() + delimiter
                        + result.getSplicingPathomechanism());
                writer.newLine();
            }
        }
    }

    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("threes-evaluator")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }


        List<SimpleResults> results = new ArrayList<>();
        // -----------------------    FOR EACH PHENOPACKET    --------------------------------------
        for (Path phenopacketPath : phenopacketPaths) {
            LOGGER.info("Reading phenopacket from '{}'", phenopacketPath);
            Phenopacket pp;
            pp = Utils.readPhenopacket(phenopacketPath);
            if (pp.getSubject().getId().isEmpty()) {
                LOGGER.warn("Phenopacket subject's ID must not be empty. Unable to continue");
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

            String entrezString = pp.getGenes(0).getId();
            if (entrezString.contains("ENTREZ:")) {
                entrezString = entrezString.substring(7);
            }
            int entrezId = Integer.parseInt(entrezString);


            // -----------------------    CREATE THE SIMULATED VCF FILE    -------------------------
            LOGGER.info("Creating simulated VCF file");
            VcfSimulator simulator = new SingleVcfSimulator(templateVcfPath);
            Path vcfPath = simulator.simulateVcfWithPhenopacket(pp);


            // -------------------------------------------------------------------------------------
            //
            // First run the analysis without splicing scores
            //
            LOGGER.info("Creating splicing-agnostic analysis");
            String sampleName = pp.getSubject().getId().replaceAll("\\s+", "_");
            InheritanceModeOptions inheritanceModeOptions = getInheritanceModeOptions(entrezId);

            List<String> phenotypesAsHpoStrings = Utils.getPresentPhenotypesAsHpoStrings(pp);

            Analysis splicingAgnosticAnalysis = exomiser.getAnalysisBuilder()
                    .genomeAssembly(GenomeAssembly.HG19)
                    .analysisMode(AnalysisMode.PASS_ONLY)
//                    .inheritanceModes(inheritanceModeOptions)
                    .vcfPath(vcfPath)
                    .probandSampleName(sampleName)
                    .frequencySources(FrequencySource.ALL_EXTERNAL_FREQ_SOURCES)
                    .pathogenicitySources(PS_NOT_SPLICING) // all the pathogenicity sources except SPLICING & TEST
                    .hpoIds(phenotypesAsHpoStrings)
                    .addFrequencyFilter(FREQ_CUTOFF)
//                    .addFrequencyFilter() // use frequencies set by inheritanceModeOptions
                    .addPathogenicityFilter(true)
                    .addInheritanceFilter()
                    .addOmimPrioritiser()
                    .addHiPhivePrioritiser()
                    .build();

            AnalysisResults splicingAgnosticResults = exomiser.run(splicingAgnosticAnalysis);

            OutputSettings agnosticSettings = OutputSettings.builder()
                    .outputFormats(EnumSet.of(OutputFormat.HTML))
                    .outputPrefix(outputPath.resolve(phenopacketPath.toFile().getName() + "_NO").toString())
                    .build();
            AnalysisResultsWriter.writeToFile(splicingAgnosticAnalysis, splicingAgnosticResults, agnosticSettings);

            //
            // Then run the analysis with splicing pathogenicity scores
            LOGGER.info("Creating splicing-aware analysis");
            Analysis splicingAwareAnalysis = exomiser.getAnalysisBuilder()
                    .genomeAssembly(GenomeAssembly.HG19)
                    .analysisMode(AnalysisMode.PASS_ONLY)
//                    .inheritanceModes(inheritanceModeOptions)
                    .vcfPath(vcfPath)
                    .probandSampleName(sampleName)
                    .frequencySources(FrequencySource.ALL_EXTERNAL_FREQ_SOURCES)
                    .pathogenicitySources(PS_W_SPLICING) // all the pathogenicity sources except TEST
                    .hpoIds(phenotypesAsHpoStrings)
                    .addFrequencyFilter(FREQ_CUTOFF)
//                    .addFrequencyFilter() // use frequencies set by inheritanceModeOptions
                    .addPathogenicityFilter(true)
                    .addInheritanceFilter()
                    .addOmimPrioritiser()
                    .addHiPhivePrioritiser()
                    .build();

            AnalysisResults splicingAwareResults = exomiser.run(splicingAwareAnalysis);

            OutputSettings awareSettings = OutputSettings.builder()
                    .outputFormats(EnumSet.of(OutputFormat.HTML))
                    .outputPrefix(outputPath.resolve(phenopacketPath.toFile().getName() + "_YES").toString())
                    .build();
            AnalysisResultsWriter.writeToFile(splicingAwareAnalysis, splicingAwareResults, awareSettings);

            // store evaluation of the analyses
            results.add(evaluateResults(pp, splicingAgnosticResults, splicingAwareResults));

        }

        writeResults(results, outputPath);

        LOGGER.info("Done");
    }

    /**
     * Get inheritance modes of diseases associated with mutations in gene with given {@code entrezGeneId}.
     *
     * @param entrezGeneId entrez id of affected gene
     * @return
     */
    private InheritanceModeOptions getInheritanceModeOptions(int entrezGeneId) {
        if (inheritanceModeAware) {
            List<Disease> diseases = priorityService.getDiseaseDataAssociatedWithGeneId(entrezGeneId);

            return InheritanceModeOptions.defaultForModes(diseases.stream()
                    .map(Disease::getInheritanceMode)
                    .map(InheritanceMode::toModeOfInheritance)
                    .flatMap(Collection::stream)
                    .distinct()
                    .toArray(ModeOfInheritance[]::new));
        } else {
            // all inheritance modes with default MAX AF
            return InheritanceModeOptions.defaults();
        }
    }

    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopackets
        if (!args.containsOption("pp")) {
            LOGGER.warn("At least one '--pp' argument for Phenopacket path must be present");
            return false;
        }
        phenopacketPaths.addAll(args.getOptionValues("pp").stream().map(Paths::get).collect(Collectors.toList()));

        //
        if (!args.containsOption("vcf")) {
            LOGGER.warn("Missing 'vcf' argument");
            return false;
        }
        templateVcfPath = Paths.get(args.getOptionValues("vcf").get(0));

        // Output file path - results
        if (!args.containsOption("output")) {
            LOGGER.warn("Missing '--output' argument");
            return false;
        }
        outputPath = Paths.get(args.getOptionValues("output").get(0));

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
