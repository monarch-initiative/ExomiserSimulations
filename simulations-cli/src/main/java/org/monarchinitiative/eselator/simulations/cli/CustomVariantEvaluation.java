package org.monarchinitiative.eselator.simulations.cli;

import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.phenopackets.schema.v1.Phenopacket;

public class CustomVariantEvaluation {

    private final VariantEvaluation variantEvaluation;

    private final String vcfAlleleInfoField;

    private final Phenopacket phenopacket;

    private final double donorScore, acceptorScore, intronScore, exonScore;

    private CustomVariantEvaluation(Builder builder) {
        this.variantEvaluation = builder.variantEvaluation;
        this.phenopacket = builder.phenopacket;
        this.donorScore = builder.donorScore;
        this.acceptorScore = builder.acceptorScore;
        this.intronScore = builder.intronScore;
        this.exonScore = builder.exonScore;
        this.vcfAlleleInfoField = builder.vcfAlleleInfoField;
    }

    public static Builder builder() {
        return new Builder();
    }

    public String getVcfAlleleInfoField() {
        return vcfAlleleInfoField;
    }

    public VariantEvaluation getVariantEvaluation() {
        return variantEvaluation;
    }

    public Phenopacket getPhenopacket() {
        return phenopacket;
    }

    public double getDonorScore() {
        return donorScore;
    }

    public double getAcceptorScore() {
        return acceptorScore;
    }

    public double getIntronScore() {
        return intronScore;
    }

    public double getExonScore() {
        return exonScore;
    }

    public static class Builder {

        private VariantEvaluation variantEvaluation;

        private Phenopacket phenopacket;

        private double donorScore, acceptorScore, intronScore, exonScore;

        private String vcfAlleleInfoField;

        private Builder() {

        }

        public Builder setVcfAlleleInfoField(String vcfAlleleInfoField) {
            this.vcfAlleleInfoField = vcfAlleleInfoField;
            return this;
        }

        public Builder setDonorScore(double donorScore) {
            this.donorScore = donorScore;
            return this;
        }

        public Builder setAcceptorScore(double acceptorScore) {
            this.acceptorScore = acceptorScore;
            return this;
        }

        public Builder setIntronScore(double intronScore) {
            this.intronScore = intronScore;
            return this;
        }

        public Builder setExonScore(double exonScore) {
            this.exonScore = exonScore;
            return this;
        }

        public Builder setVariantEvaluation(VariantEvaluation variantEvaluation) {
            this.variantEvaluation = variantEvaluation;
            return this;
        }

        public Builder setPhenopacket(Phenopacket phenopacket) {
            this.phenopacket = phenopacket;
            return this;
        }

        public CustomVariantEvaluation build() {
            return new CustomVariantEvaluation(this);
        }
    }
}

