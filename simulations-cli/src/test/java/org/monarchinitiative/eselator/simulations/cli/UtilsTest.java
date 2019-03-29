package org.monarchinitiative.eselator.simulations.cli;

import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.jupiter.api.Test;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Individual;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;

import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.CoreMatchers.*;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.monarchinitiative.eselator.simulations.cli.TestExamples.*;

class UtilsTest {




    @Test
    void functionPhenopacketToHetVariantContextWorks() {
        String jb = "JohnnyBravo";
        List<VariantContext> variantContexts = Utils.phenopacketToVariantContexts(
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
        List<VariantContext> variantContexts = Utils.phenopacketToVariantContexts(
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
        List<VariantContext> variantContexts = Utils.phenopacketToVariantContexts(
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



}