package org.monarchinitiative.threes.benchmarks.cmd.score_splicing_variants;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

class AnnotatedVariantDataTest {

    private AnnotatedVariantData data;

    @BeforeEach
    void setUp() {
        data = AnnotatedVariantData.builder()
                .rawVariantData(RawVariantData.builder()
                        .variantId("idee")
                        .txId("txidd")
                        .contig("chr4")
                        .begin(49)
                        .end(342)
                        .ref("A")
                        .alt("C")
                        .highestEffect("SYNO")
                        .effects("SYNO")
                        .cChange("234A>C")
                        .variantSource("CLN")
                        .build())
                .splicingPathogenicityData(SplicingPathogenicityData.builder()
                        .putScore("canonical_donor", .1)
                        .putScore("cryptic_donor", .2)
                        .putScore("canonical_acceptor", .3)
                        .putScore("cryptic_acceptor", .4)
                        .build())
                .build();
    }

    @Test
    void toRecord() {
        final Object[] actual = data.meltToRecord();
        Object[] expected = new Object[]{"idee", "txidd", "chr4", 49, 342, "A", "C", "SYNO", "SYNO", "234A>C", "CLN", .1, .2, .3, .4};
        assertArrayEquals(expected, actual);
    }

}