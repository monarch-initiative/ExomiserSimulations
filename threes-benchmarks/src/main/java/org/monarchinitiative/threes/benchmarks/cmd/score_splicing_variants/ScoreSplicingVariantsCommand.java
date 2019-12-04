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

        subparser.addArgument("input")
                .help("path to CSV file with one variant per line");

        subparser.addArgument("output")
                .help("where to write the output file");
    }

    private static Function<CSVRecord, RawVariantData> toRawVariantData() {
        return record -> RawVariantData.builder()
                .variantId(record.get("VARIANT_ID"))
                .txId(record.get("TX_ID"))
                .contig(record.get("CONTIG"))
                .begin(Integer.parseInt(record.get("BEGIN")))
                .end(Integer.parseInt(record.get("END")))
                .ref(record.get("REF"))
                .alt(record.get("ALT"))
                .effects(record.get("EFFECTS"))
                .cChange(record.get("C_CHANGE"))
                .variantSource(record.get("VARIANT_SOURCE"))
                .build();
    }

    @Override
    public void run(Namespace namespace) throws CommandException {
        LOGGER.info("Running `score-splicing-variants`");
        // 0 - parse CLI args
        final Path input = Paths.get(namespace.getString("input"));
        final Path output = Paths.get(namespace.getString("output"));


        // 1 - decode variants from the input file
        List<RawVariantData> variants;
        List<String> headerNames;
        try (final BufferedReader reader = Files.newBufferedReader(input);
             final CSVParser parse = CSVFormat.DEFAULT.withFirstRecordAsHeader().parse(reader)) {
            LOGGER.info("Reading variants");
            headerNames = new ArrayList<>(parse.getHeaderNames());
            // VARIANT_ID	TX_ID	CONTIG	BEGIN	END	REF	ALT	EFFECTS	C_CHANGE	VARIANT_SOURCE
            variants = StreamSupport.stream(parse.spliterator(), false)
                    .map(toRawVariantData())
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
        String[] header = headerNames.toArray(new String[0]);
        LOGGER.info("Writing results into `{}`", output);
        try (final BufferedWriter writer = Files.newBufferedWriter(output);
             final CSVPrinter csvPrinter = CSVFormat.DEFAULT.withHeader(header).print(writer)) {
            annotated.stream()
                    .map(AnnotatedVariantData::meltToRecord)
                    .forEach(values -> {
                        try {
                            csvPrinter.printRecord(values);
                        } catch (IOException e) {
                            LOGGER.warn("Error writing record `{}`", values);
                        }
                    });
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
