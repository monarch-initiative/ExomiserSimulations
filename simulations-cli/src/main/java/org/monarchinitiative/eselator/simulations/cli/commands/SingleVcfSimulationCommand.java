package org.monarchinitiative.eselator.simulations.cli.commands;

import org.monarchinitiative.eselator.simulations.cli.Utils;
import org.monarchinitiative.eselator.simulations.cli.simulators.SingleVcfSimulator;
import org.monarchinitiative.eselator.simulations.cli.simulators.VcfSimulator;
import org.monarchinitiative.exomiser.core.Exomiser;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisMode;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.writers.AnalysisResultsWriter;
import org.monarchinitiative.exomiser.core.writers.OutputFormat;
import org.monarchinitiative.exomiser.core.writers.OutputSettings;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Phenotype;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@Component
public class SingleVcfSimulationCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(SingleVcfSimulationCommand.class);

    private static final float MAX_FREQ = 0.001f;

    private final Exomiser exomiser;

    private Path templateVcfPath;

    private Path phenopacketFilePath;

    private Path outputPrefixPath;

    public SingleVcfSimulationCommand(Exomiser exomiser) {
        this.exomiser = exomiser;
    }


    /**
     * Only present (non-negated) {@link Phenotype}s are reported
     *
     * @param pp {@link Phenopacket} describing the proband
     * @return list of HPO id strings representing subject's phenotype
     */
    static List<String> getPresentPhenotypesAsHpoStrings(Phenopacket pp) {
        return pp.getPhenotypesList().stream()
                .filter(p -> !p.getAbsent())
                .map(p -> p.getType().getId())
                .collect(Collectors.toList());
    }

    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("single-vcf-simulation")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }

        // -----------------------    READ PHENOPACKET    ------------------------------------------
        LOGGER.info("Reading phenopacket from '{}'", phenopacketFilePath);
        Phenopacket pp;
        pp = Utils.readPhenopacket(phenopacketFilePath);
        if (pp.getSubject().getId().isEmpty()) {
            LOGGER.warn("Phenopacket subject's ID must not be empty. Unable to continue");
            System.exit(1);
        }

        // -----------------------    CREATE THE SIMULATED VCF FILE    -----------------------------
        LOGGER.info("Creating simulated VCF file");
        VcfSimulator simulator = new SingleVcfSimulator(templateVcfPath);
        Path vcfPath = simulator.simulateVcfWithPhenopacket(pp);


        // -----------------------    FORGE EXOMISER ANALYSIS    -----------------------------------
        Set<FrequencySource> frequencySources = new HashSet<>(FrequencySource.FREQUENCY_SOURCE_MAP.values());
        LOGGER.info("Creating analysis");
        Analysis analysis = exomiser.getAnalysisBuilder()
                .analysisMode(AnalysisMode.PASS_ONLY)
                .vcfPath(vcfPath)
                .probandSampleName(pp.getSubject().getId())
                .genomeAssembly(GenomeAssembly.HG19)
                .frequencySources(frequencySources)
                .addFrequencyFilter(MAX_FREQ)
                .hpoIds(getPresentPhenotypesAsHpoStrings(pp))
                .addOmimPrioritiser()
                .addHiPhivePrioritiser()
                .build();

        // -----------------------    RUN THE ANALYSIS AND WRITE THE RESULTS    --------------------
        LOGGER.info("Running the analysis");
        final AnalysisResults results = exomiser.run(analysis);

        OutputSettings settings = OutputSettings.builder()
                .outputFormats(EnumSet.allOf(OutputFormat.class))
                .outputPrefix(outputPrefixPath.toString())
                .build();
        AnalysisResultsWriter.writeToFile(analysis, results, settings);

        LOGGER.info("Done!");
    }

    /**
     * Parse the command line arguments and return false if anything is missing. Errors are also logged.
     */
    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopacket
        if (!args.containsOption("pp")) {
            LOGGER.warn("Missing '--pp' argument for Phenopacket");
            return false;
        }
        phenopacketFilePath = Paths.get(args.getOptionValues("pp").get(0));

        // Template VCF file
        if (!args.containsOption("vcf")) {
            LOGGER.warn("Missing '--vcf' argument for template VCF file path");
            return false;
        }
        templateVcfPath = Paths.get(args.getOptionValues("vcf").get(0));

        // Output path prefix
        if (!args.containsOption("output")) {
            LOGGER.warn("Missing '--output' argument for where to write the results");
            return false;
        }
        outputPrefixPath = Paths.get(args.getOptionValues("output").get(0));

        return true;
    }
}
