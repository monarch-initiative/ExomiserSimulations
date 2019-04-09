package org.monarchinitiative.eselator.simulations.cli;

import org.monarchinitiative.exomiser.autoconfigure.EnableExomiser;
import org.monarchinitiative.exomiser.core.genome.dao.splicing.SplicingParameters;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.annotation.Bean;
import org.springframework.core.env.Environment;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Objects;

@EnableExomiser
@SpringBootApplication
public class Play {

    public static void main(String[] args) {
        SpringApplication.run(Play.class, args);
    }

    @Bean
    public Path hexamerFilePath(Environment environment) {
        return Paths.get(Objects.requireNonNull(environment.getProperty("hexamer.scores.path")));
    }
}
