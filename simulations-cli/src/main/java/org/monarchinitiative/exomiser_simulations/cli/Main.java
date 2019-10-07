package org.monarchinitiative.exomiser_simulations.cli;

import org.monarchinitiative.exomiser.autoconfigure.EnableExomiser;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

//@EnableExomiser
@SpringBootApplication
public class Main {

    public static void main(String[] args) {
        try {
            SpringApplication.run(Main.class, args);
        } catch (Throwable t) {
            t.printStackTrace(); // this will show when the user has entered the wrong path for input files etc.
        }
    }
}
