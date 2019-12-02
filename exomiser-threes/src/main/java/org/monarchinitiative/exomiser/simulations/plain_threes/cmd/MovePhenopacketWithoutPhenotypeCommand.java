package org.monarchinitiative.exomiser.simulations.plain_threes.cmd;

import org.monarchinitiative.exomiser.simulations.plain_threes.Utils;
import org.phenopackets.schema.v1.Phenopacket;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * This commands finds phenopackets that contain no phenotype term (HPO id) and which would thus crash Exomiser analysis.
 * <p>
 * These phenopackets without phenotype terms are moved into the {@code noPhenotypeDirName} directory that is created in
 * the {@code phenopacketDir} directory.
 */
@Component
public class MovePhenopacketWithoutPhenotypeCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(MovePhenopacketWithoutPhenotypeCommand.class);

    private String noPhenotypeDirName = "_no_phenotype";

    private File phenopacketDir;

    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("move-phenopackets-without-phenotype")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }

        final File[] pps = Objects.requireNonNull(phenopacketDir.listFiles(f -> f.getName().endsWith(".json")));

        Arrays.sort(pps, Comparator.comparing(File::getName));
        List<File> withoutPhenotype = new ArrayList<>();

        for (File ppfile : pps) {
            final Phenopacket phenopacket = Utils.readPhenopacket(ppfile.toPath());
            System.out.println(String.format("%s\t%d terms", ppfile.getName(), phenopacket.getPhenotypicFeaturesCount()));
            if (phenopacket.getPhenotypicFeaturesList().isEmpty()) {
                withoutPhenotype.add(ppfile);
            }
        }

        System.out.println("\n# ----------------------     EMPTY     ----------------------\n");
        System.out.println(String.format("%d phenopackets without phenotype", withoutPhenotype.size()));

        Path noPhenotype = phenopacketDir.toPath().resolve(noPhenotypeDirName);
        Files.createDirectories(noPhenotype); // ensure the directory exists
        for (File file : withoutPhenotype) {
            // move
            Files.move(file.toPath(), noPhenotype.resolve(file.getName()));
        }
    }

    private boolean parseCliArgs(ApplicationArguments args) {

        if (!args.containsOption("pp-dir")) {
            LOGGER.warn("Missing 'pp-dir' argument");
            return false;
        }
        phenopacketDir = new File(args.getOptionValues("pp-dir").get(0));

        if (!phenopacketDir.isDirectory()) {
            LOGGER.warn("Argument 'pp-dir' does not point to directory");
            return false;
        }
        return true;
    }
}
