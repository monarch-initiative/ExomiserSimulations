package org.monarchinitiative.exomiser.simulations.plain_threes;

import com.google.common.collect.ImmutableMap;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.phenopackets.schema.v1.Phenopacket;

import java.util.HashMap;
import java.util.Map;

public class CustomVariantEvaluation {

    private final VariantEvaluation variantEvaluation;

    private final String vcfAlleleInfoField;

    private final Phenopacket phenopacket;

    private final ImmutableMap<String, Double> scoreMap;

    private CustomVariantEvaluation(Builder builder) {
        this.variantEvaluation = builder.variantEvaluation;
        this.phenopacket = builder.phenopacket;
        this.vcfAlleleInfoField = builder.vcfAlleleInfoField;
        this.scoreMap = ImmutableMap.copyOf(builder.scoreMap);
    }

    public static Builder builder() {
        return new Builder();
    }

    public ImmutableMap<String, Double> getScoreMap() {
        return scoreMap;
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


    public static class Builder {

        private VariantEvaluation variantEvaluation;

        private Phenopacket phenopacket;

        private String vcfAlleleInfoField;

        private Map<String, Double> scoreMap = new HashMap<>();

        private Builder() {

        }

        public Builder setVcfAlleleInfoField(String vcfAlleleInfoField) {
            this.vcfAlleleInfoField = vcfAlleleInfoField;
            return this;
        }

        public Builder putScore(String key, double value) {
            this.scoreMap.put(key, value);
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

