package org.monarchinitiative.exomiser.simulations.cli.commands;

import org.monarchinitiative.exomiser.core.Exomiser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 *
 */
@Component
public class ThreesEvaluatorCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(ThreesEvaluatorCommand.class);

    /**
     * List of paths pointing to phenopackets that should be evaluated.
     */
    private final List<Path> phenopacketPaths = new ArrayList<>();

    private final Exomiser exomiser;

    /**
     * Path to directory where output will be directed.
     */
    private Path outputPath;

    public ThreesEvaluatorCommand(Exomiser exomiser) {
        this.exomiser = exomiser;
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

        //
        LOGGER.info("Evaluation is running");

        //
        LOGGER.info("Done");
    }

    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopackets
        if (!args.containsOption("pp")) {
            LOGGER.warn("At least one '--pp' argument for Phenopacket path must be present");
            return false;
        }
        phenopacketPaths.addAll(args.getOptionValues("pp").stream().map(Paths::get).collect(Collectors.toList()));

        // Output file path - results
        if (!args.containsOption("output")) {
            LOGGER.warn("Missing '--output' argument");
            return false;
        }
        outputPath = Paths.get(args.getOptionValues("output").get(0));

        return true;
    }
}
