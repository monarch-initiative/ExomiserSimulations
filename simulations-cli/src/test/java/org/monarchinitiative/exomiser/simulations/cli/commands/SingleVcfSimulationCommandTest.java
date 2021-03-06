package org.monarchinitiative.exomiser.simulations.cli.commands;

import org.junit.jupiter.api.Test;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.PhenotypicFeature;

import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.CoreMatchers.*;
import static org.monarchinitiative.exomiser.simulations.cli.TestExamples.ontologyClass;

class SingleVcfSimulationCommandTest {


    @Test
    void phenotypesFromPhenopacketAreExtractedToListOfHpoIds() {
        final Phenopacket pp = Phenopacket.newBuilder()
                .addPhenotypicFeatures(PhenotypicFeature.newBuilder().setType(ontologyClass("HP:0000252", "Microcephaly")).build())
                .addPhenotypicFeatures(PhenotypicFeature.newBuilder().setType(ontologyClass("HP:0001382", "Joint hypermobility")).build())
                .addPhenotypicFeatures(PhenotypicFeature.newBuilder().setType(ontologyClass("HP:0003357", "Thymic hormone decreased")).setNegated(true).build())
                .build();

        final List<String> hpos = SingleVcfSimulationCommand.getPresentPhenotypesAsHpoStrings(pp);

        assertThat(hpos.size(), is(2));
        assertThat(hpos, hasItems("HP:0000252", "HP:0001382"));
    }
}