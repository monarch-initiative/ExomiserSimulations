package org.monarchinitiative.threes.benchmarks;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import net.sourceforge.argparse4j.inf.Subparsers;
import org.monarchinitiative.threes.autoconfigure.EnableThrees;
import org.monarchinitiative.threes.benchmarks.cmd.ClinvarScorerCommand;
import org.monarchinitiative.threes.benchmarks.cmd.Command;
import org.monarchinitiative.threes.benchmarks.cmd.ScorePhenopacketsCommand;
import org.monarchinitiative.threes.benchmarks.cmd.score_splicing_variants.ScoreSplicingVariantsCommand;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.builder.SpringApplicationBuilder;
import org.springframework.context.ConfigurableApplicationContext;

import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Properties;

@EnableThrees
@SpringBootApplication
public class Main {
    private static final Logger LOGGER = LoggerFactory.getLogger(Main.class);

    private static final String EPILOG = "       くまくま━━━━━━ヽ（ ・(ｪ)・ ）ノ━━━━━━ !!!";

    public static void main(String[] args) throws Exception {

        // 1. define CLI interface
        ArgumentParser parser = ArgumentParsers.newFor("java -jar threes-benchmarks.jar").build();
        parser.description("Benchmarks of 3S code");

        // - we require threes configuration/properties
        parser.addArgument("-c", "--config")
                .required(true)
                .metavar("/path/to/application.properties")
                .help("path to Spring configuration file");

        Subparsers subparsers = parser.addSubparsers();
        ClinvarScorerCommand.setupSubparsers(subparsers);
        ScorePhenopacketsCommand.setupSubparsers(subparsers);
        ScoreSplicingVariantsCommand.setupSubparsers(subparsers);

        parser.defaultHelp(true);
        parser.epilog(EPILOG);

        Namespace namespace = null;
        try {
            namespace = parser.parseArgs(args);
        } catch (ArgumentParserException e) {
            parser.handleError(e);
            System.exit(1);
        }

        Properties properties = new Properties();
        try (InputStream is = Files.newInputStream(Paths.get(namespace.getString("config")))) {
            properties.load(is);
        }

        //  2. bootstrap the app
        try (ConfigurableApplicationContext appContext = new SpringApplicationBuilder(Main.class)
                .properties(properties)
                .run()) {

            // 3. get the selected command and run it
            Command command;
            String cmdName = namespace.get("cmd");
            switch (cmdName) {
                case "clinvar-scorer":
                    command = appContext.getBean(ClinvarScorerCommand.class);
                    break;
                case "score-phenopackets":
                    command = appContext.getBean(ScorePhenopacketsCommand.class);
                    break;
                case "score-splicing-variants":
                    command = appContext.getBean(ScoreSplicingVariantsCommand.class);
                    break;
                default:
                    LOGGER.warn("Unknown command '{}'", cmdName);
                    System.exit(1);
                    return; // unreachable, but still required
            }

            command.run(namespace);
            LOGGER.info("Done!");
        }
    }
}
