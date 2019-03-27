package org.monarchinitiative.eselator.simulations;

import org.monarchinitiative.exomiser.autoconfigure.EnableExomiser;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

@SpringBootApplication
@EnableExomiser
public class Play {

    public static void main(String[] args) {
        SpringApplication.run(Play.class, args);
    }
}
