package org.monarchinitiative.eselator.simulations.app;

import org.monarchinitiative.exomiser.core.Exomiser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

@Component
public class SimulationsRunner implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(SimulationsRunner.class);

    private final Exomiser exomiser;

    public SimulationsRunner(Exomiser exomiser) {
        this.exomiser = exomiser;
    }

    public void run(ApplicationArguments args) throws Exception {
        LOGGER.info("###########################################################################");
        LOGGER.info("#                                                                         #");
        LOGGER.info("#                     EXOMISER SET UP SUCCESSFULLY!!                      #");
        LOGGER.info("#                                                                         #");
        LOGGER.info("###########################################################################");
    }
}
