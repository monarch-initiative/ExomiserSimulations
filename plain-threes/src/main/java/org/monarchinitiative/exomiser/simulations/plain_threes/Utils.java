package org.monarchinitiative.exomiser.simulations.plain_threes;

import com.google.protobuf.util.JsonFormat;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.monarchinitiative.threes.core.model.SplicingTranscript;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.OntologyClass;
import org.phenopackets.schema.v1.core.PhenotypicFeature;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public final class Utils {

    public static final OntologyClass HET = OntologyClass.newBuilder().setId("GENO:0000135").setLabel("heterozygous").build();

    public static final OntologyClass HOM_ALT = OntologyClass.newBuilder().setId("GENO:0000136").setLabel("homozygous").build();

    public static final OntologyClass HEMIZYGOUS = OntologyClass.newBuilder().setId("GENO:0000134").setLabel("hemizygous").build();

    private static final Logger LOGGER = LoggerFactory.getLogger(Utils.class);

    private static final JsonFormat.Parser PARSER = JsonFormat.parser();

    private static final Pattern INFO_BIFIELD = Pattern.compile("(\\w+)=(-?[\\S]+)");

    private Utils() {
        // private no-op
    }

    /**
     * @param phenopacketFilePath {@link Path} to Phenopacket with data in JSON format
     * @return decoded {@link Phenopacket}
     * @throws IOException in case of <code>phenopacketFilePath</code> I/O errors, if JSON data has invalid format, and
     *                     if
     */
    public static Phenopacket readPhenopacket(Path phenopacketFilePath) throws IOException {
        Phenopacket.Builder builder = Phenopacket.newBuilder();
        try (BufferedReader reader = Files.newBufferedReader(phenopacketFilePath)) {
            String json = reader.lines().collect(Collectors.joining());
            PARSER.merge(json, builder);
        }
        return builder.build();
    }

    /**
     * Map {@link Variant}s to {@link VariantContext}s. Genotypes in variant contexts are modified so that they
     * will contain the subject's id.
     */
    public static List<VariantContext> phenopacketToVariantContexts(List<Variant> variants, String subjectId) {
        List<VariantContext> variantContexts = new ArrayList<>();
        for (Variant variant : variants) {

            switch (variant.getAlleleCase()) {
                case ALLELE_NOT_SET:
                case SPDI_ALLELE:
                case ISCN_ALLELE:
                case HGVS_ALLELE:
                default:
                    LOGGER.warn("Variant data are not stored in VCF format, but as {}", variant.getAlleleCase());
                    continue;
                case VCF_ALLELE:
                    // continue execution
            }

            final VcfAllele vcfAllele = variant.getVcfAllele();

            // here the ref allele is always at 0, alt is at idx 1
            List<Allele> allAlleles = new ArrayList<>(2);
            allAlleles.add(Allele.create(vcfAllele.getRef(), true));
            allAlleles.add(Allele.create(vcfAllele.getAlt()));

            OntologyClass zygosity = variant.getZygosity();
            GenotypeBuilder genotypeBuilder = new GenotypeBuilder()
                    .name(subjectId);

            if (zygosity.equals(HET)) {
                // 1x REF + 1x ALT
                genotypeBuilder.alleles(Arrays.asList(allAlleles.get(0), allAlleles.get(1)));
            } else if (zygosity.equals(HOM_ALT)) {
                // 2x ALT
                genotypeBuilder.alleles(Arrays.asList(allAlleles.get(1), allAlleles.get(1)));
            } else if (zygosity.equals(HEMIZYGOUS)) {
                genotypeBuilder.alleles(Collections.singletonList(allAlleles.get(1)));
            } else {
                LOGGER.warn("Unknown genotype '{}'. Tried HET, HOM_ALT, HEMIZYGOUS", zygosity);
                continue;
            }

            // INFO fields

            Map<String, Object> infoFields = new HashMap<>();
            for (String s : vcfAllele.getInfo()
                    .replaceAll("\\s+", "_")
                    .split(";")) {
                Matcher bifield = INFO_BIFIELD.matcher(s);
                if (bifield.matches()) {
                    infoFields.put(bifield.group(1), bifield.group(2));
                } else {
                    infoFields.put(s, null);
                }
            }

            VariantContext vc = new VariantContextBuilder()
                    // we are working with hg19 usually. Contigs are prepended with 'chr' there
                    .chr(vcfAllele.getChr().startsWith("chr") ? vcfAllele.getChr() : "chr" + vcfAllele.getChr())
                    .start(vcfAllele.getPos())
                    .computeEndFromAlleles(allAlleles, vcfAllele.getPos())
                    .alleles(allAlleles)
                    .genotypes(genotypeBuilder.make())
                    .attributes(infoFields)
                    .noID()
                    .make();

            variantContexts.add(vc);
        }
        return variantContexts;
    }


    /**
     * Only present (non-negated) {@link PhenotypicFeature}s are reported
     *
     * @param pp {@link Phenopacket} describing the proband
     * @return list of HPO id strings representing subject's phenotype
     */
    public static List<String> getPresentPhenotypesAsHpoStrings(Phenopacket pp) {
        return pp.getPhenotypicFeaturesList().stream()
                .filter(p -> !p.getNegated())
                .map(p -> p.getType().getId())
                .collect(Collectors.toList());
    }


    public static String getSplicingPathomechanism(List<Variant> variants) {
        if (variants.isEmpty()) {
            return "None";
        }
        List<String> pathomechanisms = new ArrayList<>(3);
        for (Variant variant : variants) {
            if (variant.getAlleleCase() == Variant.AlleleCase.VCF_ALLELE) {
                final VcfAllele vcfAllele = variant.getVcfAllele();
                // alleleInfo might look like:
                // `VCLASS=splicing;PATHOMECHANISM=splicing|3ss|disrupted;CONSEQUENCE=Alternative/cryptic 3' splice site`
                final String alleleInfo = vcfAllele.getInfo();
                final String[] tokens = alleleInfo.split(";");
                for (String token : tokens) {
                    final String[] fields = token.split("=");
                    if (fields[0].equals("PATHOMECHANISM")) {
                        pathomechanisms.add(fields[1]);
                    }
                }
            }
        }

        return String.join(";", pathomechanisms);
    }
    /**
     * @return {@link Comparator} where the highest priority transcript is the one with the smallest integer value of the
     * central part of the accession ID. The central part of the <em>NM_004004.2</em> is <em>004004</em>.
     */
    public static Comparator<? super SplicingTranscript> transcriptPriorityComparator() {
        return (l, r) -> {
            int leftInt = getCentralInt(l.getAccessionId());
            int rightInt = getCentralInt(r.getAccessionId());
            return Integer.compare(leftInt, rightInt);
        };
    }

    /**
     * @param accId String like <em>NM_004004.2</em>
     * @return <em>4004</em> - the central integer part of the transcript accession id string
     */
    private static int getCentralInt(String accId) {
        if (accId.matches("NM_\\d+.?\\d*")) {
            int dotIdx = accId.indexOf(".");
            String central = accId.substring(3, dotIdx);
            return Integer.parseInt(central);
        }

        throw new RuntimeException("Weird accession ID " + accId);
    }

}
