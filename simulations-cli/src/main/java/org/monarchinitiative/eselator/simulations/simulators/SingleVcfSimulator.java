package org.monarchinitiative.eselator.simulations.simulators;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.OntologyClass;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

/**
 * Simulator that injects variants defined from {@link Phenopacket} among variants present in single VCF file.
 */
public class SingleVcfSimulator implements VcfSimulator {

    public static final OntologyClass HET = OntologyClass.newBuilder().setId("GENO:0000135").setLabel("heterozygous").build();

    public static final OntologyClass HOM_ALT = OntologyClass.newBuilder().setId("GENO:0000136").setLabel("homozygous").build();

    public static final OntologyClass HEMIZYGOUS = OntologyClass.newBuilder().setId("GENO:0000134").setLabel("hemizygous").build();

    private static final Pattern INFO_BIFIELD = Pattern.compile("(\\w+)=(-?[\\w.]+)");


    private static final Logger LOGGER = LoggerFactory.getLogger(SingleVcfSimulator.class);

    private final Path templateVcfPath;

    /**
     * @param templateVcfPath {@link Path} to possibly un-indexed VCF file
     */
    public SingleVcfSimulator(Path templateVcfPath) {
        this.templateVcfPath = templateVcfPath;
    }


    static List<VariantContext> phenopacketToVariantContexts(Phenopacket phenopacket, String sampleId) {
        List<VariantContext> variants = new ArrayList<>();
        for (Variant variant : phenopacket.getVariantsList()) {

            VcfAllele vcfAllele = variant.getVcfAllele();
            if (!vcfAllele.isInitialized()) { // VcfAllele is oneof filed, thus not always present
                LOGGER.warn("VcfAllele is not present in variant {}", variant);
                continue;
            }

            // here the ref allele is always at 0, alt is at idx 1
            List<Allele> allAlleles = new ArrayList<>(2);
            allAlleles.add(Allele.create(vcfAllele.getRef(), true));
            allAlleles.add(Allele.create(vcfAllele.getAlt()));

            OntologyClass genotype = variant.getGenotype();
            GenotypeBuilder genotypeBuilder = new GenotypeBuilder()
                    .name(sampleId);

            if (genotype.equals(HET)) {
                // 1x REF + 1x ALT
                genotypeBuilder.alleles(Arrays.asList(allAlleles.get(0), allAlleles.get(1)));
            } else if (genotype.equals(HOM_ALT)) {
                // 2x ALT
                genotypeBuilder.alleles(Arrays.asList(allAlleles.get(1), allAlleles.get(1)));
            } else if (genotype.equals(HEMIZYGOUS)) {
                genotypeBuilder.alleles(Collections.singletonList(allAlleles.get(1)));
            } else {
                LOGGER.warn("Unknown genotype '{}'. Tried HET, HOM_ALT, HEMIZYGOUS", genotype);
                continue;
            }

            // INFO fields

            Map<String, Object> infoFields = new HashMap<>();
            for (String s : vcfAllele.getInfo().split(";")) {
                Matcher bifield = INFO_BIFIELD.matcher(s);
                if (bifield.matches()) {
                    infoFields.put(bifield.group(1), bifield.group(2));
                } else {
                    infoFields.put(s, null);
                }
            }


            VariantContext vc = new VariantContextBuilder()
                    .chr("chr" + vcfAllele.getChr())
                    .start(vcfAllele.getPos())
                    .computeEndFromAlleles(allAlleles, vcfAllele.getPos())
                    .alleles(allAlleles)
                    .genotypes(genotypeBuilder.make())
                    .attributes(infoFields)
                    .noID()
                    .make();

            variants.add(vc);
        }
        return variants;
    }

    @Override
    public Path simulateVcfWithPhenopacket(Phenopacket phenopacket) throws IOException {
        File outPath = File.createTempFile("single-vcf-simulator-", ".vcf");
        outPath.deleteOnExit();

        VCFFileReader reader = new VCFFileReader(templateVcfPath, false);
        VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(outPath)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();
        try (reader; writer) {
            LOGGER.info("Reading file {}", templateVcfPath);
            VCFHeader fileHeader = reader.getFileHeader();
            writer.writeHeader(fileHeader);

            ArrayList<String> samples = fileHeader.getSampleNamesInOrder();
            List<VariantContext> injected = phenopacketToVariantContexts(phenopacket, samples.get(0));

           Stream.concat(reader.iterator().stream(), injected.stream())
                    .sorted(new VariantContextComparator(fileHeader.getContigLines()))
                    .forEach(writer::add);
        }
        return outPath.toPath();
    }
}
