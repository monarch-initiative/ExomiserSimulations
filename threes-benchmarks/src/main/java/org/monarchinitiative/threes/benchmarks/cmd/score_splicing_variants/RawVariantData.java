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
    private final String symbol;
    private final int donorOffset;
    private final int acceptorOffset;

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
        symbol = builder.symbol;
        offset = builder.offset;
        donorOffset = builder.donorOffset;
        acceptorOffset = builder.acceptorOffset;
    }

    public static Builder builder() {
        return new Builder();
    }

    public int getDonorOffset() {
        return donorOffset;
    }

    public int getAcceptorOffset() {
        return acceptorOffset;
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

    public String getSymbol() {
        return symbol;
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
                donorOffset == that.donorOffset &&
                acceptorOffset == that.acceptorOffset &&
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
                Objects.equals(clz, that.clz) &&
                Objects.equals(symbol, that.symbol);
    }

    @Override
    public int hashCode() {
        return Objects.hash(variantId, txId, contig, begin, end, ref, alt, highestEffect, effects, cChange, variantSource, culprit, clz, offset, symbol, donorOffset, acceptorOffset);
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
                ", symbol='" + symbol + '\'' +
                ", donorOffset=" + donorOffset +
                ", acceptorOffset=" + acceptorOffset +
                '}';
    }

    public static final class Builder {
        private String variantId;
        private String txId;
        private String symbol;
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
        private int donorOffset;
        private int acceptorOffset;

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

        public Builder symbol(String symbol) {
            this.symbol = symbol;
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

        public Builder donorOffset(int offset) {
            this.donorOffset = offset;
            return this;
        }

        public Builder acceptorOffset(int offset) {
            this.acceptorOffset = offset;
            return this;
        }

        public RawVariantData build() {
            return new RawVariantData(this);
        }
    }
}
