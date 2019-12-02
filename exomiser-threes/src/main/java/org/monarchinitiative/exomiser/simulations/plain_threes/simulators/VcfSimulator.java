package org.monarchinitiative.exomiser.simulations.plain_threes.simulators;

import org.phenopackets.schema.v1.Phenopacket;

import java.io.IOException;
import java.nio.file.Path;

public interface VcfSimulator {


    Path simulateVcfWithPhenopacket(Phenopacket phenopacket) throws IOException;


}
