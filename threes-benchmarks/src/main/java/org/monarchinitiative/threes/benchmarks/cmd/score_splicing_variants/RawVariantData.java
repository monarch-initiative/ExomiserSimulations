package org.monarchinitiative.threes.benchmarks.cmd.score_splicing_variants;

import java.util.Objects;

class RawVariantData {

    private final String variantId;
    private final String txId;
    private final String contig;
    private final int begin;
    private final int end;
    private final String ref;
    private final String alt;
    private final String highestEffect;
    private final String effects;
    private final String cChange;
    private final String variantSource;
    private final String culprit;
    private final String clz;
    private final int offset;

    private RawVariantData(Builder builder) {
        variantId = builder.variantId;
        txId = builder.txId;
        contig = builder.contig;
        begin = builder.begin;
        end = builder.end;
        ref = builder.ref;
        alt = builder.alt;
        highestEffect = builder.highestEffect;
        effects = builder.effects;
        cChange = builder.cChange;
        variantSource = builder.variantSource;
        culprit = builder.culprit;
        clz = builder.clz;
        offset = builder.offset;
    }

    public static Builder builder() {
        return new Builder();
    }

    public String getCulprit() {
        return culprit;
    }

    public String getClz() {
        return clz;
    }

    public int getOffset() {
        return offset;
    }

    public String getHighestEffect() {
        return highestEffect;
    }

    public String getVariantId() {
        return variantId;
    }

    public String getTxId() {
        return txId;
    }

    public String getContig() {
        return contig;
    }

    public int getBegin() {
        return begin;
    }

    public int getEnd() {
        return end;
    }

    public String getRef() {
        return ref;
    }

    public String getAlt() {
        return alt;
    }

    public String getEffects() {
        return effects;
    }

    public String getcChange() {
        return cChange;
    }

    public String getVariantSource() {
        return variantSource;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        RawVariantData that = (RawVariantData) o;
        return begin == that.begin &&
                end == that.end &&
                offset == that.offset &&
                Objects.equals(variantId, that.variantId) &&
                Objects.equals(txId, that.txId) &&
                Objects.equals(contig, that.contig) &&
                Objects.equals(ref, that.ref) &&
                Objects.equals(alt, that.alt) &&
                Objects.equals(highestEffect, that.highestEffect) &&
                Objects.equals(effects, that.effects) &&
                Objects.equals(cChange, that.cChange) &&
                Objects.equals(variantSource, that.variantSource) &&
                Objects.equals(culprit, that.culprit) &&
                Objects.equals(clz, that.clz);
    }

    @Override
    public int hashCode() {
        return Objects.hash(variantId, txId, contig, begin, end, ref, alt, highestEffect, effects, cChange, variantSource, culprit, clz, offset);
    }

    @Override
    public String toString() {
        return "RawVariantData{" +
                "variantId='" + variantId + '\'' +
                ", txId='" + txId + '\'' +
                ", contig='" + contig + '\'' +
                ", begin=" + begin +
                ", end=" + end +
                ", ref='" + ref + '\'' +
                ", alt='" + alt + '\'' +
                ", highestEffect='" + highestEffect + '\'' +
                ", effects='" + effects + '\'' +
                ", cChange='" + cChange + '\'' +
                ", variantSource='" + variantSource + '\'' +
                ", culprit='" + culprit + '\'' +
                ", clz='" + clz + '\'' +
                ", offset=" + offset +
                '}';
    }

    public static final class Builder {
        private String variantId;
        private String txId;
        private String contig;
        private int begin;
        private int end;
        private String ref;
        private String alt;
        private String highestEffect;
        private String effects;
        private String cChange;
        private String variantSource;
        private String culprit;
        private String clz;
        private int offset;

        private Builder() {
        }

        public Builder variantId(String variantId) {
            this.variantId = variantId;
            return this;
        }

        public Builder txId(String txId) {
            this.txId = txId;
            return this;
        }

        public Builder contig(String contig) {
            this.contig = contig;
            return this;
        }

        public Builder begin(int begin) {
            this.begin = begin;
            return this;
        }

        public Builder end(int end) {
            this.end = end;
            return this;
        }

        public Builder ref(String ref) {
            this.ref = ref;
            return this;
        }

        public Builder alt(String alt) {
            this.alt = alt;
            return this;
        }

        public Builder highestEffect(String highestEffect) {
            this.highestEffect = highestEffect;
            return this;
        }

        public Builder effects(String effects) {
            this.effects = effects;
            return this;
        }

        public Builder cChange(String cChange) {
            this.cChange = cChange;
            return this;
        }

        public Builder variantSource(String variantSource) {
            this.variantSource = variantSource;
            return this;
        }

        public Builder culprit(String culprit) {
            this.culprit = culprit;
            return this;
        }

        public Builder clz(String clz) {
            this.clz = clz;
            return this;
        }

        public Builder offset(int offset) {
            this.offset = offset;
            return this;
        }

        public RawVariantData build() {
            return new RawVariantData(this);
        }
    }
}
