package org.monarchinitiative.eselator.simulations.cli;

import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Individual;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;

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
import static org.monarchinitiative.eselator.simulations.cli.TestExamples.*;

class SingleVcfSimulatorTest {

    private static final Path TEST_VCF_PATH = Paths.get(SingleVcfSimulatorTest.class.getResource("GIAB_NIST7035_3vars.vcf").getFile());

    private SingleVcfSimulator instance;

    /**
     * Make Phenopacket for sample <code>NIST7035</code> - <code>chr1:787400C>T HET AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.08</code>
     *
     * @return blah
     */
    private static Phenopacket makePhenopacketWithHetVariant(Variant variant, Individual subject) {
        return Phenopacket.newBuilder()
                .setSubject(subject)
                .addVariants(variant)
                .setMetaData(getMetaData(getGenotypeOntologyResource()))
                .build();
    }

    private static Individual individual(String individualId) {
        return Individual.newBuilder()
                .setId(individualId)
                .build();
    }

    private static Variant hetVariant() {
        return Variant.newBuilder()
                .setVcfAllele(VcfAllele.newBuilder()
                        .setChr("chr1")
                        .setPos(787400)
                        .setRef("C")
                        .setAlt("T")
                        .setInfo("AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.08")
                        .build())
                .setGenotype(HET)
                .build();
    }

    private static Variant homAltVariant() {
        return Variant.newBuilder()
                .setVcfAllele(VcfAllele.newBuilder()
                        .setChr("chr1")
                        .setPos(787400)
                        .setRef("C")
                        .setAlt("T")
                        .setInfo("AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.08")
                        .build())
                .setGenotype(HOM_ALT)
                .build();
    }

    private static Variant hemizygousVariant() {
        return Variant.newBuilder()
                .setVcfAllele(VcfAllele.newBuilder()
                        .setChr("chr1")
                        .setPos(787400)
                        .setRef("C")
                        .setAlt("T")
                        .setInfo("AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.08")
                        .build())
                .setGenotype(HEMIZYGOUS)
                .build();
    }

    @BeforeEach
    void setUp() throws Exception {
        instance = new SingleVcfSimulator(TEST_VCF_PATH);
    }

    @Test
    void initializationWorks() throws IOException {
        Phenopacket packet = makePhenopacketWithHetVariant(hetVariant(), individual("ID"));
        VariantContext expected = SingleVcfSimulator.phenopacketToVariantContexts(packet).get(0);

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

    @Test
    void functionPhenopacketToHetVariantContextWorks() {
        String jb = "JohnnyBravo";
        List<VariantContext> variantContexts = SingleVcfSimulator.phenopacketToVariantContexts(
                makePhenopacketWithHetVariant(hetVariant(), individual(jb)));

        assertThat(variantContexts.size(), is(1));
        VariantContext vc = variantContexts.get(0);

        assertThat(vc.getContig(), is("chr1"));
        assertThat(vc.getStart(), is(787400));
        assertThat(vc.getAlleles().size(), is(2));
        assertThat(vc.getReference().getBaseString(), is("C"));
        assertThat(vc.getAlternateAllele(0).getBaseString(), is("T"));
        assertTrue(vc.hasGenotype(jb));
        assertThat(vc.getGenotype(jb).getType(), is(GenotypeType.HET));
    }


    @Test
    void functionPhenopacketToHomaltVariantContextWorks() {
        String jb = "JohnnyBravo";
        List<VariantContext> variantContexts = SingleVcfSimulator.phenopacketToVariantContexts(
                makePhenopacketWithHetVariant(homAltVariant(), individual(jb)));

        assertThat(variantContexts.size(), is(1));
        VariantContext vc = variantContexts.get(0);

        assertThat(vc.getContig(), is("chr1"));
        assertThat(vc.getStart(), is(787400));
        assertThat(vc.getAlleles().size(), is(2));
        assertThat(vc.getReference().getBaseString(), is("C"));
        assertThat(vc.getAlternateAllele(0).getBaseString(), is("T"));
        assertTrue(vc.hasGenotype(jb));
        assertThat(vc.getGenotype(jb).getType(), is(GenotypeType.HOM_VAR));
    }

    @Test
    void functionPhenopacketToHemizygousVariantContextWorks() {
        String jb = "JohnnyBravo";
        List<VariantContext> variantContexts = SingleVcfSimulator.phenopacketToVariantContexts(
                makePhenopacketWithHetVariant(hemizygousVariant(), individual(jb)));

        assertThat(variantContexts.size(), is(1));
        VariantContext vc = variantContexts.get(0);

        assertThat(vc.getContig(), is("chr1"));
        assertThat(vc.getStart(), is(787400));
        assertThat(vc.getAlleles().size(), is(2));
        assertThat(vc.getReference().getBaseString(), is("C"));
        assertThat(vc.getAlternateAllele(0).getBaseString(), is("T"));
        assertTrue(vc.hasGenotype(jb));
        assertThat(vc.getGenotype(jb).getType(), is(GenotypeType.HOM_VAR));
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