package org.monarchinitiative.exomiser.simulations.cli;

import org.monarchinitiative.exomiser.autoconfigure.EnableExomiser;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

@EnableExomiser
@SpringBootApplication
public class Play {

    public static void main(String[] args) {
        SpringApplication.run(Play.class, args);
    }

}
