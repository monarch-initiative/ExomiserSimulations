package org.monarchinitiative.eselator.simulations.simulators;

import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.junit.Before;
import org.junit.Test;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Individual;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static org.monarchinitiative.eselator.simulations.TestUtils.*;

public class SingleVcfSimulatorTest {

    private static final Path TEST_VCF_PATH = Paths.get(SingleVcfSimulatorTest.class.getResource("GIAB_NIST7035_3vars.vcf").getFile());

    private static final String SAMPLE_ID = "NIST7035";

    private SingleVcfSimulator instance;

    /**
     * Make Phenopacket for sample <code>NIST7035</code> - <code>chr1:787400C>T HET AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=32.08</code>
     *
     * @return blah
     */
    private static Phenopacket makePhenopacketWithHetVariant(Variant variant) {
        return Phenopacket.newBuilder()
                .addVariants(variant)
                .setMetaData(getMetaData(getGenotypeOntologyResource()))
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

    @Before
    public void setUp() throws Exception {
        instance = new SingleVcfSimulator(TEST_VCF_PATH);
    }

    @Test
    public void initializationWorks() throws IOException {
        Phenopacket packet = makePhenopacketWithHetVariant(hetVariant());
        VariantContext expected = SingleVcfSimulator.phenopacketToVariantContexts(packet, SAMPLE_ID).get(0);

        Path path = instance.simulateVcfWithPhenopacket(packet);
        List<VariantContext> variants = new ArrayList<>();
        try (VCFFileReader reader = new VCFFileReader(path.toFile(), false)) {
            variants.addAll(reader.iterator().toList());
        }
        // TODO - continue here
        assertThat(variants.size(), is(4));
        assertTrue(variants.stream().anyMatch(vc -> vc.getContig().equals(expected.getContig())
                && vc.getStart() == 787400
                && vc.getEnd() == 787400
                && vc.getReference().getBaseString().equals("C")
                && vc.getAlternateAllele(0).getBaseString().equals("T")));
    }

    @Test
    public void functionPhenopacketToHetVariantContextWorks() {
        List<VariantContext> variantContexts = SingleVcfSimulator.phenopacketToVariantContexts(makePhenopacketWithHetVariant(hetVariant()), SAMPLE_ID);

        assertThat(variantContexts.size(), is(1));
        VariantContext vc = variantContexts.get(0);

        assertThat(vc.getContig(), is("chr1"));
        assertThat(vc.getStart(), is(787400));
        assertThat(vc.getAlleles().size(), is(2));
        assertThat(vc.getReference().getBaseString(), is("C"));
        assertThat(vc.getAlternateAllele(0).getBaseString(), is("T"));
        assertTrue(vc.hasGenotype("NIST7035"));
        assertThat(vc.getGenotype("NIST7035").getType(), is(GenotypeType.HET));
    }


    @Test
    public void functionPhenopacketToHomaltVariantContextWorks() {
        List<VariantContext> variantContexts = SingleVcfSimulator.phenopacketToVariantContexts(makePhenopacketWithHetVariant(homAltVariant()), SAMPLE_ID);

        assertThat(variantContexts.size(), is(1));
        VariantContext vc = variantContexts.get(0);

        assertThat(vc.getContig(), is("chr1"));
        assertThat(vc.getStart(), is(787400));
        assertThat(vc.getAlleles().size(), is(2));
        assertThat(vc.getReference().getBaseString(), is("C"));
        assertThat(vc.getAlternateAllele(0).getBaseString(), is("T"));
        assertTrue(vc.hasGenotype("NIST7035"));
        assertThat(vc.getGenotype("NIST7035").getType(), is(GenotypeType.HOM_VAR));
    }

    @Test
    public void functionPhenopacketToHemizygousVariantContextWorks() {
        List<VariantContext> variantContexts = SingleVcfSimulator.phenopacketToVariantContexts(makePhenopacketWithHetVariant(hemizygousVariant()), SAMPLE_ID);

        assertThat(variantContexts.size(), is(1));
        VariantContext vc = variantContexts.get(0);

        assertThat(vc.getContig(), is("chr1"));
        assertThat(vc.getStart(), is(787400));
        assertThat(vc.getAlleles().size(), is(2));
        assertThat(vc.getReference().getBaseString(), is("C"));
        assertThat(vc.getAlternateAllele(0).getBaseString(), is("T"));
        assertTrue(vc.hasGenotype("NIST7035"));
        assertThat(vc.getGenotype("NIST7035").getType(), is(GenotypeType.HOM_VAR));
    }
}