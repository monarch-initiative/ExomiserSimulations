package org.monarchinitiative.threes.benchmarks.cmd.score_splicing_variants;

import de.charite.compbio.jannovar.reference.*;
import net.sourceforge.argparse4j.inf.Namespace;
import net.sourceforge.argparse4j.inf.Subparser;
import net.sourceforge.argparse4j.inf.Subparsers;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.csv.CSVRecord;
import org.monarchinitiative.threes.benchmarks.cmd.Command;
import org.monarchinitiative.threes.benchmarks.cmd.CommandException;
import org.monarchinitiative.threes.core.data.SplicingTranscriptSource;
import org.monarchinitiative.threes.core.model.SplicingTranscript;
import org.monarchinitiative.threes.core.scoring.SplicingAnnotator;
import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;
import xyz.ielis.hyperutil.reference.fasta.GenomeSequenceAccessor;
import xyz.ielis.hyperutil.reference.fasta.SequenceInterval;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

@Component
public class ScoreSplicingVariantsCommand extends Command {

    private static final Logger LOGGER = LoggerFactory.getLogger(ScoreSplicingVariantsCommand.class);
    private static final int TX_PADDING = 150;
    private final GenomeSequenceAccessor genomeSequenceAccessor;
    private final SplicingTranscriptSource splicingTranscriptSource;
    private final SplicingAnnotator splicingAnnotator;

    public ScoreSplicingVariantsCommand(GenomeSequenceAccessor genomeSequenceAccessor,
                                        SplicingTranscriptSource splicingTranscriptSource,
                                        SplicingAnnotator splicingAnnotator) {
        this.genomeSequenceAccessor = genomeSequenceAccessor;
        this.splicingTranscriptSource = splicingTranscriptSource;
        this.splicingAnnotator = splicingAnnotator;
    }


    public static void setupSubparsers(Subparsers subparsers) {
        Subparser subparser = subparsers.addParser("score-splicing-variants")
                .setDefault("cmd", "score-splicing-variants")
                .help("score splicing variants from CSV file");

        subparser.addArgument("-t", "--input-type")
                .type(String.class)
                .choices(List.of("default", "overall"))
                .setDefault("default")
                .help("we have 2 input types which have different columns, which one are we getting now?");

        subparser.addArgument("input")
                .help("path to CSV file with one variant per line");

        subparser.addArgument("output")
                .help("where to write the output file");
    }

    /**
     * Get decoder function for given `inputType`.
     *
     * @param inputType String with input type
     * @return {@link Function} for decoding CSV line into raw variant data
     */
    private static Function<CSVRecord, RawVariantData> parseInputType(String inputType) {
        switch (inputType) {
            case "default":
                return record -> RawVariantData.builder()
                        .variantId(record.get("VARIANT_ID"))
                        .txId(record.get("TX_ID"))
                        .contig(record.get("CONTIG"))
                        .begin(Integer.parseInt(record.get("BEGIN")))
                        .end(Integer.parseInt(record.get("END")))
                        .ref(record.get("REF"))
                        .alt(record.get("ALT"))
                        .highestEffect(record.get("HIGHEST_EFFECT"))
                        .effects(record.get("EFFECTS"))
                        .cChange(record.get("C_CHANGE"))
                        .variantSource(record.get("VARIANT_SOURCE"))
                        .build();
            case "overall":
                return record -> RawVariantData.builder()
                        .variantId(record.get("variant_id"))
                        .contig(record.get("contig"))
                        .begin(Integer.parseInt(record.get("begin")))
                        .end(Integer.parseInt(record.get("end")))
                        .ref(record.get("ref"))
                        .alt(record.get("alt"))
                        .txId(record.get("tx_id"))
                        .culprit(record.get("culprit"))
                        .variantSource(record.get("source"))
                        .clz(record.get("clz"))
                        .offset(Integer.parseInt(record.get("offset")))
                        .symbol(record.get("symbol"))
                        .donorOffset(Integer.parseInt(record.get("donor_offset")))
                        .acceptorOffset(Integer.parseInt(record.get("acceptor_offset")))
                        .build();
            default:
                LOGGER.warn("Unknown input type `{}`", inputType);
                throw new RuntimeException(String.format("Unknown input type `%s`", inputType));
        }
    }

    private static Consumer<AnnotatedVariantData> writeOutVariant(CSVPrinter printer, String inputType) {
        switch (inputType) {
            case "default":
                return data -> {
                    // here the header looks like:
                    // variant_id,tx_id,contig,begin,end,
                    // ref,alt,highest_effect,effects,
                    // c_change,source,
                    // canonical_donor,cryptic_donor,
                    // canonical_acceptor,cryptic_acceptor
                    try {
                        final RawVariantData rvd = data.getRawVariantData();
                        final SplicingPathogenicityData spd = data.getSplicingPathogenicityData();
                        printer.printRecord(rvd.getVariantId(), rvd.getTxId(), rvd.getContig(), rvd.getBegin(), rvd.getEnd(),
                                rvd.getRef(), rvd.getAlt(), rvd.getHighestEffect(), String.join(",", rvd.getEffects()),
                                rvd.getcChange(), rvd.getVariantSource(),
                                spd.getOrDefault("canonical_donor", Double.NaN), spd.getOrDefault("cryptic_donor", Double.NaN),
                                spd.getOrDefault("canonical_acceptor", Double.NaN), spd.getOrDefault("cryptic_acceptor", Double.NaN));
                    } catch (IOException e) {
                        LOGGER.warn("Error writing record `{}`", data);
                    }
                };
            case "overall":
                // here the header looks like:
                // variant_id,contig,begin,end,
                // ref,alt,tx_id,
                // culprit,source,clz,offset,
                // symbol,donor_offset,acceptor_offset,
                // canonical_donor,cryptic_donor,
                // canonical_acceptor,cryptic_acceptor
                return data -> {
                    try {
                        final RawVariantData rvd = data.getRawVariantData();
                        final SplicingPathogenicityData spd = data.getSplicingPathogenicityData();
                        printer.printRecord(rvd.getVariantId(), rvd.getContig(), rvd.getBegin(), rvd.getEnd(),
                                rvd.getRef(), rvd.getAlt(), rvd.getTxId(),
                                rvd.getCulprit(), rvd.getVariantSource(), rvd.getClz(), rvd.getOffset(),
                                rvd.getSymbol(), rvd.getDonorOffset(), rvd.getAcceptorOffset(),
                                spd.getOrDefault("canonical_donor", Double.NaN), spd.getOrDefault("cryptic_donor", Double.NaN),
                                spd.getOrDefault("canonical_acceptor", Double.NaN), spd.getOrDefault("cryptic_acceptor", Double.NaN));
                    } catch (IOException e) {
                        LOGGER.warn("Error writing record `{}`", data);
                    }
                };
            default:
                LOGGER.warn("Unknown input type `{}`", inputType);
                throw new RuntimeException(String.format("Unknown input type `%s`", inputType));
        }
    }

    @Override
    public void run(Namespace namespace) throws CommandException {
        LOGGER.info("Running `score-splicing-variants`");
        // 0 - parse CLI args
        final Path input = Paths.get(namespace.getString("input"));
        final Path output = Paths.get(namespace.getString("output"));
        final String inputType = namespace.getString("input_type");

        // 1 - decode variants from the input file
        List<RawVariantData> variants;
        List<String> headerNames;
        Function<CSVRecord, RawVariantData> decoder = parseInputType(inputType);

        try (final BufferedReader reader = Files.newBufferedReader(input);
             final CSVParser csvParser = CSVFormat.DEFAULT
                     .withFirstRecordAsHeader()
                     .parse(reader)) {
            LOGGER.info("Reading variants");
            // header names depend on `input_type`
            headerNames = new ArrayList<>(csvParser.getHeaderNames());
            variants = StreamSupport.stream(csvParser.spliterator(), false)
                    .map(decoder)
                    .collect(Collectors.toList());
            LOGGER.info("Read {} variants", variants.size());
        } catch (IOException e) {
            throw new CommandException(e);
        }

        // 2 - group variants by transcripts and perform annotation
        final Map<String, List<RawVariantData>> variantsByTx = variants.stream().collect(Collectors.groupingBy(RawVariantData::getTxId));
        LOGGER.info("Annotating variants on {} transcripts", variantsByTx.size());
        final List<AnnotatedVariantData> annotated = variantsByTx.entrySet().parallelStream()
                .map(annotateVariantsOfTranscript())
                .flatMap(Collection::stream)
                .collect(Collectors.toList());

        // 3 - write the annotated variants into a CSV file
        headerNames.addAll(AnnotatedVariantData.SCORER_NAMES);
        String[] header = headerNames.stream().map(String::toLowerCase).toArray(String[]::new);
        LOGGER.info("Writing {} results into `{}`", annotated.size(), output);
        try (final BufferedWriter writer = Files.newBufferedWriter(output);
             final CSVPrinter csvPrinter = CSVFormat.DEFAULT.withHeader(header).print(writer)) {
            annotated.forEach(writeOutVariant(csvPrinter, inputType));
        } catch (IOException e) {
            throw new CommandException(e);
        }
    }

    /**
     * Annotate variants of a transcript. The entry key is transcript accession id, the value is list of variants.
     *
     * @return function for annotation as described above
     */
    private Function<Map.Entry<String, List<RawVariantData>>, Collection<AnnotatedVariantData>> annotateVariantsOfTranscript() {
        return entry -> {
            // fetch transcript
            final Optional<SplicingTranscript> sto = splicingTranscriptSource.fetchTranscriptByAccession(entry.getKey(), genomeSequenceAccessor.getReferenceDictionary());
            if (sto.isEmpty()) {
                LOGGER.warn("Splicing transcript for `{}` not found", entry.getKey());
                return Collections.emptySet();
            }
            final SplicingTranscript st = sto.get();

            // fetch reference sequence
            final GenomeInterval txRegion = st.getTxRegionCoordinates().withMorePadding(TX_PADDING, TX_PADDING);
            final Optional<SequenceInterval> sio = genomeSequenceAccessor.fetchSequence(txRegion);
            if (sio.isEmpty()) {
                LOGGER.warn("Unable to fetch sequence for `{}` from coordinates `{}`", entry.getKey(), txRegion);
                return Collections.emptySet();
            }
            final SequenceInterval si = sio.get();

            return entry.getValue().stream()
                    .map(variant -> {
                        final GenomePosition position = new GenomePosition(genomeSequenceAccessor.getReferenceDictionary(), Strand.FWD,
                                genomeSequenceAccessor.getReferenceDictionary().getContigNameToID().get(variant.getContig()),
                                variant.getBegin(),
                                PositionType.ZERO_BASED);
                        final GenomeVariant gv = new GenomeVariant(position, variant.getRef(), variant.getAlt());
                        final SplicingPathogenicityData spd = splicingAnnotator.evaluate(gv, st, si);

                        return AnnotatedVariantData.builder()
                                .rawVariantData(variant)
                                .splicingPathogenicityData(spd)
                                .build();
                    })
                    .collect(Collectors.toSet());
        };
    }

}
