package org.monarchinitiative.exomiser.simulations.cli.simulators;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.monarchinitiative.exomiser.simulations.cli.TestExamples;
import org.monarchinitiative.exomiser.simulations.cli.Utils;
import org.phenopackets.schema.v1.Phenopacket;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static org.hamcrest.CoreMatchers.hasItem;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.monarchinitiative.exomiser.simulations.cli.TestExamples.hetVariant;
import static org.monarchinitiative.exomiser.simulations.cli.TestExamples.individual;

class SingleVcfSimulatorTest {

    private static final Path TEST_VCF_PATH = Paths.get(SingleVcfSimulatorTest.class.getResource("GIAB_NIST7035_3vars.vcf").getFile());

    private SingleVcfSimulator instance;


    @BeforeEach
    void setUp() throws Exception {
        instance = new SingleVcfSimulator(TEST_VCF_PATH);
    }

    @Test
    void initializationWorks() throws IOException {
        Phenopacket packet = TestExamples.makePhenopacketWithHetVariant(hetVariant(), individual("ID"));
        VariantContext expected = Utils.phenopacketToVariantContexts(packet.getVariantsList(), "ID").get(0);

        Path path = instance.simulateVcfWithPhenopacket(packet);
        List<VariantContext> variants = new ArrayList<>();
        try (VCFFileReader reader = new VCFFileReader(path.toFile(), false)) {
            variants.addAll(reader.iterator().toList());
        }

        assertThat(variants.size(), is(4));
        assertTrue(variants.stream().anyMatch(vc -> vc.getContig().equals(expected.getContig())
                && vc.getStart() == 787400
                && vc.getEnd() == 787400
                && vc.getReference().getBaseString().equals("C")
                && vc.getAlternateAllele(0).getBaseString().equals("T")));
    }


    @ParameterizedTest
    @ValueSource(strings = {"Johnny", "Donna"})
    void changingSampleIdInVcfHeaderWorks(String sampleName) {
        final VCFHeader header;
        try (VCFFileReader reader = new VCFFileReader(
                new File(getClass().getResource("GIAB_NIST7035.vcf").getFile()), false)) {
            header = reader.getFileHeader();
        }
        final ArrayList<String> names = header.getSampleNamesInOrder();
        assertThat(names, hasItem("NIST7035"));


        final VCFHeader vcfHeader = SingleVcfSimulator.updateHeaderWithPhenopacketSample(header, sampleName);

        final ArrayList<String> updatedNames = vcfHeader.getSampleNamesInOrder();
        assertThat(updatedNames.size(), is(1));
        assertThat(updatedNames, hasItem(sampleName));
    }
}