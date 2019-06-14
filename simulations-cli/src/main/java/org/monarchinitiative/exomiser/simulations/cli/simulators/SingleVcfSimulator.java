package org.monarchinitiative.exomiser.simulations.cli.simulators;

import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.monarchinitiative.exomiser.simulations.cli.Utils;
import org.phenopackets.schema.v1.Phenopacket;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Simulator that injects variants defined from {@link Phenopacket} among variants present in single VCF file.
 */
public class SingleVcfSimulator implements VcfSimulator {


    private static final Logger LOGGER = LoggerFactory.getLogger(SingleVcfSimulator.class);

    private final Path templateVcfPath;

    /**
     * @param templateVcfPath {@link Path} to possibly un-indexed VCF file
     */
    public SingleVcfSimulator(Path templateVcfPath) {
        this.templateVcfPath = templateVcfPath;
    }


    static VCFHeader updateHeaderWithPhenopacketSample(VCFHeader original, String sampleId) {
        return new VCFHeader(original.getMetaDataInSortedOrder(), Collections.singleton(sampleId));
    }

    static UnaryOperator<VariantContext> changeSampleNameInGenotypes(final String sampleId) {
        return vc -> {
            final VariantContextBuilder vcb = new VariantContextBuilder(vc)
                    .noGenotypes() // remove present genotypes and then add updated
                    .genotypes(vc.getGenotypes().stream()
                            .map(gt -> new GenotypeBuilder(gt).name(sampleId).make()) // change sample Id on individual genotypes
                            .collect(Collectors.toList()));

            return vcb.make();
        };
    }

    @Override
    public Path simulateVcfWithPhenopacket(Phenopacket phenopacket) throws IOException {

        final String sampleId = phenopacket.getSubject().getId();
        // we create a temporary VCF file for Exomiser analysis
        final File outPath = File.createTempFile("single-vcf-simulators-" + sampleId + "-", ".vcf");
        outPath.deleteOnExit();

        try (VCFFileReader reader = new VCFFileReader(templateVcfPath, false);
             VariantContextWriter writer = new VariantContextWriterBuilder()
                     .setOutputFile(outPath)
                     .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                     .unsetOption(Options.INDEX_ON_THE_FLY)
                     .build()) {
            LOGGER.info("Reading file {}", templateVcfPath);
            VCFHeader fileHeader = reader.getFileHeader();
            fileHeader = updateHeaderWithPhenopacketSample(fileHeader, phenopacket.getSubject().getId());
            writer.writeHeader(fileHeader);

            List<VariantContext> injected = Utils.phenopacketToVariantContexts(phenopacket);

            AtomicInteger cnt = new AtomicInteger();
            Stream.concat(reader.iterator().stream(), injected.stream())
                    .map(changeSampleNameInGenotypes(sampleId))
                    .sorted(new VariantContextComparator(fileHeader.getContigLines()))
                    .peek(vc -> cnt.incrementAndGet())
                    .forEach(writer::add);
            LOGGER.info("Created VCF containing {} variants", cnt.get());
        }
        return outPath.toPath();
    }
}
