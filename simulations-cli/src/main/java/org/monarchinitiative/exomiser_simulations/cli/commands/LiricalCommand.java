package org.monarchinitiative.exomiser_simulations.cli.commands;

import org.monarchinitiative.exomiser_simulations.cli.Utils;
import org.monarchinitiative.exomiser_simulations.cli.simulators.SingleVcfSimulator;
import org.monarchinitiative.exomiser_simulations.cli.simulators.VcfSimulator;
import org.monarchinitiative.exomiser.core.Exomiser;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisMode;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Phenotype;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@Component
public class LiricalCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(LiricalCommand.class);

    /**
     * HARDCODED FOR NOW
     */
    private static final float MAX_FREQ = 0.001f;

    private final Exomiser exomiser;

    private Path templateVcfPath;

    private Path phenopacketDirectoryPath;

    public LiricalCommand(Exomiser exomiser) {
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
                .filter(p -> !p.getNegated())
                .map(p -> p.getType().getId())
                .collect(Collectors.toList());
    }

    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("lirical")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }

        // -----------------------    FOR EACH PHENOPACKET    --------------------------------------
        File[] fileArray = phenopacketDirectoryPath.toFile().listFiles();
        if (fileArray == null) {
            LOGGER.warn("Phenopacket file array is null. This should not happen");
            return;
        }
        List<File> phenopackets = Arrays.asList(fileArray);
        for (File phenopacketFilePath : phenopackets) {
            // -----------------------    READ PHENOPACKET    --------------------------------------
            LOGGER.info("Reading phenopacket from '{}'", phenopacketFilePath);
            Phenopacket pp;
            pp = Utils.readPhenopacket(phenopacketFilePath.toPath());
            if (pp.getSubject().getId().isEmpty()) {
                LOGGER.warn("Phenopacket subject's ID must not be empty. Unable to continue");
                System.exit(1);
            }

            // -----------------------    CREATE THE SIMULATED VCF FILE    -------------------------
            LOGGER.info("Creating simulated VCF file");
            VcfSimulator simulator = new SingleVcfSimulator(templateVcfPath);
            Path vcfPath = simulator.simulateVcfWithPhenopacket(pp);


            // -----------------------    FORGE EXOMISER ANALYSIS    -------------------------------
            Set<FrequencySource> frequencySources = new HashSet<>(FrequencySource.ALL_EXTERNAL_FREQ_SOURCES);
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


            // -----------------------    RUN THE ANALYSIS AND WRITE THE RESULTS    ----------------
            LOGGER.info("Running the analysis");
            final AnalysisResults results = exomiser.run(analysis);
            // TODO(pnrobinson) - evaluate
            System.out.println(results);
        }

    }

    /**
     * Parse the command line arguments and return false if anything is missing. Errors are also logged.
     */
    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopacket
        if (!args.containsOption("pp")) {
            LOGGER.warn("Missing '--pp' argument for Phenopacket directory path");
            return false;
        }
        phenopacketDirectoryPath = Paths.get(args.getOptionValues("pp").get(0));
        if (!phenopacketDirectoryPath.toFile().isDirectory()) {
            LOGGER.warn("Path '{}' does not point to existing directory", phenopacketDirectoryPath);
            return false;
        }

        // Template VCF file
        if (!args.containsOption("vcf")) {
            LOGGER.warn("Missing '--vcf' argument for template VCF file path");
            return false;
        }
        templateVcfPath = Paths.get(args.getOptionValues("vcf").get(0));

        return true;
    }
}
