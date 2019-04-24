package org.monarchinitiative.exomiser_simulations.cli.simulators;

import org.phenopackets.schema.v1.Phenopacket;

import java.io.IOException;
import java.nio.file.Path;

public interface VcfSimulator {


    Path simulateVcfWithPhenopacket(Phenopacket phenoPacket) throws IOException;


}
