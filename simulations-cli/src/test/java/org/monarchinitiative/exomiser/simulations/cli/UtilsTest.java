package org.monarchinitiative.exomiser.simulations.cli;

import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.jupiter.api.Test;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;

import java.util.ArrayList;
import java.util.List;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.jupiter.api.Assertions.assertTrue;

class UtilsTest {




    @Test
    void functionPhenopacketToHetVariantContextWorks() {
        String jb = "JohnnyBravo";
        final Phenopacket pp = TestExamples.makePhenopacketWithHetVariant(TestExamples.hetVariant(), TestExamples.individual(jb));
        List<VariantContext> variantContexts = Utils.phenopacketToVariantContexts(
                pp.getVariantsList(), jb);

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
        final Phenopacket pp = TestExamples.makePhenopacketWithHetVariant(TestExamples.homAltVariant(), TestExamples.individual(jb));
        List<VariantContext> variantContexts = Utils.phenopacketToVariantContexts(
                pp.getVariantsList(), jb);

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
        final Phenopacket pp = TestExamples.makePhenopacketWithHetVariant(TestExamples.hemizygousVariant(), TestExamples.individual(jb));
        List<VariantContext> variantContexts = Utils.phenopacketToVariantContexts(
                pp.getVariantsList(), jb);

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
    void getSplicingPathomechanismFromEmptyVariantList() {
        List<Variant> variants = new ArrayList<>();

        String pathomechanism = Utils.getSplicingPathomechanism(variants);
        assertThat(pathomechanism, is("None"));
    }

    @Test
    void getSplicingPathomechanismWorksWithSingleVariant() {
        List<Variant> variants = new ArrayList<>();

        variants.add(Variant.newBuilder()
                .setVcfAllele(VcfAllele.newBuilder()
                        .setGenomeAssembly("GRCh37").setChr("chr1").setPos(123456).setRef("C").setAlt("G")
                        .setInfo("VCLASS=splicing;PATHOMECHANISM=splicing|3ss|disrupted;CONSEQUENCE=Alternative/cryptic 3' splice site")
                        .build())
                .build());

        String pathomechanism = Utils.getSplicingPathomechanism(variants);
        assertThat(pathomechanism, is("splicing|3ss|disrupted"));
    }

    @Test
    void getSplicingPathomechanismWorksWithMultipleVariants() {
        List<Variant> variants = new ArrayList<>();

        variants.add(Variant.newBuilder()
                .setVcfAllele(VcfAllele.newBuilder()
                        .setGenomeAssembly("GRCh37").setChr("chr1").setPos(123456).setRef("C").setAlt("G")
                        .setInfo("VCLASS=splicing;PATHOMECHANISM=splicing|3ss|disrupted;CONSEQUENCE=Alternative/cryptic 3' splice site")
                        .build())
                .build());
        variants.add(Variant.newBuilder()
                .setVcfAllele(VcfAllele.newBuilder()
                        .setGenomeAssembly("GRCh37").setChr("chr1").setPos(123456).setRef("C").setAlt("G")
                        .setInfo("VCLASS=splicing;PATHOMECHANISM=splicing|5css|activated;CONSEQUENCE=Alternative/cryptic 3' splice site")
                        .build())
                .build());

        String pathomechanism = Utils.getSplicingPathomechanism(variants);
        assertThat(pathomechanism, is("splicing|3ss|disrupted;splicing|5css|activated"));
    }
}