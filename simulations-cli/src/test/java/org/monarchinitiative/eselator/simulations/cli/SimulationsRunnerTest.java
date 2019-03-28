package org.monarchinitiative.eselator.simulations.cli;

import org.junit.jupiter.api.Test;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Phenotype;

import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.CoreMatchers.*;
import static org.monarchinitiative.eselator.simulations.cli.TestExamples.ontologyClass;

class SimulationsRunnerTest {


    @Test
    void phenotypesFromPhenopacketAreExtractedToListOfHpoIds() {
        final Phenopacket pp = Phenopacket.newBuilder()
                .addPhenotypes(Phenotype.newBuilder().setType(ontologyClass("HP:0000252", "Microcephaly")).build())
                .addPhenotypes(Phenotype.newBuilder().setType(ontologyClass("HP:0001382", "Joint hypermobility")).build())
                .addPhenotypes(Phenotype.newBuilder().setType(ontologyClass("HP:0003357", "Thymic hormone decreased")).setNegated(true).build())
                .build();

        final List<String> hpos = SimulationsRunner.getPresentPhenotypesAsHpoStrings(pp);

        assertThat(hpos.size(), is(2));
        assertThat(hpos, hasItems("HP:0000252", "HP:0001382"));
    }
}