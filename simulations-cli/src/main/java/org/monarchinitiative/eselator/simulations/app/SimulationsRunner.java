package org.monarchinitiative.eselator.simulations.app;

import org.monarchinitiative.eselator.simulations.simulators.SingleVcfSimulator;
import org.monarchinitiative.eselator.simulations.simulators.VcfSimulator;
import org.monarchinitiative.exomiser.core.Exomiser;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.io.PhenopacketFormat;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

@Component
public class SimulationsRunner implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(SimulationsRunner.class);

    private final Exomiser exomiser;

    private Path templateVcfPath;

    private Path phenopacketFilePath;

    public SimulationsRunner(Exomiser exomiser) {
        this.exomiser = exomiser;
    }

    public void run(ApplicationArguments args) throws Exception {
        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }

        LOGGER.info("Reading phenopacket from '{}'", phenopacketFilePath);
        Phenopacket pp;
        try (BufferedReader reader = Files.newBufferedReader(phenopacketFilePath)) {
            String json = reader.lines().collect(Collectors.joining());
            pp = PhenopacketFormat.fromJson(json);
        }

        LOGGER.info("Creating simulated VCF file");
        VcfSimulator simulator = new SingleVcfSimulator(templateVcfPath);
        Path vcfPath = simulator.simulateVcfWithPhenopacket(pp);

        LOGGER.info("Creating analysis");
        Analysis analysis = exomiser.getAnalysisBuilder()
                .vcfPath(vcfPath)
//                TODO - add the analysis steps
                .build();

        LOGGER.info("Running the analysis");

        final AnalysisResults results = exomiser.run(analysis);
    }

    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopacket
        if (!args.containsOption("pp")) {
            LOGGER.warn("Missing '--pp' argument for Phenopacket");
            return false;
        }
        List<String> pps = args.getOptionValues("pp");
        phenopacketFilePath = Paths.get(pps.get(0));

        // Template VCF file
        if (!args.containsOption("vcf")) {
            LOGGER.warn("Missing '--vcf' argument for template VCF file path");
            return false;
        }
        List<String> vcfs = args.getOptionValues("vcf");
        templateVcfPath = Paths.get(vcfs.get(0));

        return true;
    }
}
