package org.monarchinitiative.exomiser.simulations.cli.simulators;

import org.phenopackets.schema.v1.Phenopacket;

import java.io.IOException;
import java.nio.file.Path;

/**
 * This simulator creates a VCF file based on multiple VCF files. It will be determined how exactly to do this..
 */
public class MultiVcfSimulator implements VcfSimulator{

    @Override
    public Path simulateVcfWithPhenopacket(Phenopacket phenopacket) throws IOException {
        // TODO - implement
        return null;
    }
}
