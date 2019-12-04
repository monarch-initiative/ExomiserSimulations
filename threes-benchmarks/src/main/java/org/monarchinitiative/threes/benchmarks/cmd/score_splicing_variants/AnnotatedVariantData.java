package org.monarchinitiative.threes.benchmarks.cmd.score_splicing_variants;

import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;

import java.util.List;
import java.util.Objects;

class AnnotatedVariantData {

    private static final List<String> SCORER_NAMES = List.of("canonical_donor", "cryptic_donor", "canonical_acceptor", "cryptic_acceptor");

    private final RawVariantData rawVariantData;
    private final SplicingPathogenicityData splicingPathogenicityData;

    private AnnotatedVariantData(Builder builder) {
        rawVariantData = builder.rawVariantData;
        splicingPathogenicityData = builder.splicingPathogenicityData;
    }

    static Builder builder() {
        return new Builder();
    }

    public RawVariantData getRawVariantData() {
        return rawVariantData;
    }

    public SplicingPathogenicityData getSplicingPathogenicityData() {
        return splicingPathogenicityData;
    }

    Object[] meltToRecord() {
        Object[] record = new Object[14];
        record[0] = rawVariantData.getVariantId();
        record[1] = rawVariantData.getTxId();
        record[2] = rawVariantData.getContig();
        record[3] = rawVariantData.getBegin();
        record[4] = rawVariantData.getEnd();
        record[5] = rawVariantData.getRef();
        record[6] = rawVariantData.getAlt();
        record[7] = rawVariantData.getEffects();
        record[8] = rawVariantData.getcChange();
        record[9] = rawVariantData.getVariantSource();
        record[10] = splicingPathogenicityData.getOrDefault(SCORER_NAMES.get(0), Double.NaN);
        record[11] = splicingPathogenicityData.getOrDefault(SCORER_NAMES.get(1), Double.NaN);
        record[12] = splicingPathogenicityData.getOrDefault(SCORER_NAMES.get(2), Double.NaN);
        record[13] = splicingPathogenicityData.getOrDefault(SCORER_NAMES.get(3), Double.NaN);

        return record;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        AnnotatedVariantData that = (AnnotatedVariantData) o;
        return Objects.equals(rawVariantData, that.rawVariantData) &&
                Objects.equals(splicingPathogenicityData, that.splicingPathogenicityData);
    }

    @Override
    public int hashCode() {
        return Objects.hash(rawVariantData, splicingPathogenicityData);
    }

    @Override
    public String toString() {
        return "AnnotatedVariantData{" +
                "rawVariantData=" + rawVariantData +
                ", splicingPathogenicityData=" + splicingPathogenicityData +
                '}';
    }

    public static final class Builder {
        private RawVariantData rawVariantData;
        private SplicingPathogenicityData splicingPathogenicityData;

        private Builder() {
        }

        public Builder rawVariantData(RawVariantData rawVariantData) {
            this.rawVariantData = rawVariantData;
            return this;
        }

        public Builder splicingPathogenicityData(SplicingPathogenicityData splicingPathogenicityData) {
            this.splicingPathogenicityData = splicingPathogenicityData;
            return this;
        }

        public AnnotatedVariantData build() {
            return new AnnotatedVariantData(this);
        }
    }
}
