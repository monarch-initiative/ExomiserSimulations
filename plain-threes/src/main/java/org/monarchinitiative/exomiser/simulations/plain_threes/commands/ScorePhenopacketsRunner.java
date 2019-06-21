package org.monarchinitiative.exomiser.simulations.plain_threes.commands;

import org.monarchinitiative.exomiser.simulations.plain_threes.Utils;
import org.monarchinitiative.threes.core.data.SplicingTranscriptSource;
import org.monarchinitiative.threes.core.model.GenomeCoordinates;
import org.monarchinitiative.threes.core.model.SequenceInterval;
import org.monarchinitiative.threes.core.model.SplicingTranscript;
import org.monarchinitiative.threes.core.model.SplicingVariant;
import org.monarchinitiative.threes.core.reference.fasta.GenomeSequenceAccessor;
import org.monarchinitiative.threes.core.scoring.ScoringStrategy;
import org.monarchinitiative.threes.core.scoring.SplicingEvaluator;
import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * This runner implements command `--score-phenopackets`.
 * <p>
 * The runner takes directory with phenopackets (`--pp-dir=...`) as well as paths to individual phenopackets
 * (multiple `--pp=...`) as input. Variants from phenopackets are evaluated by 3S splicing evaluator and results are
 * written into TSV file (`--output-scores`).
 * <p>
 * <b>These properties have to be specified:</b>
 * <ul>
 * <li><code>`--threes.data-directory`</code> - path to directory with 3S databases & genome FASTA file</li>
 * <li><code>`--threes.genome-assembly`</code> - genome assembly - choose from {hg19, hg38}</li>
 * <li><code>`--threes.data-version`</code> - exomiser-like data version</li>
 * <li><code>`--threes.transcript-source`</code> - Jannovar transcript source - choose from {ucsc, refseq, ensembl}
 * </ul>
 * The properties are usually defined in <em>application.properties</em> file
 * </p>
 * <p>
 * In addition, you have to provide above mentioned arguments:
 * <ul>
 * <li><code>`--pp-dir`</code> - path to directory with JSON files corresponding to Phenopackets</li>
 * <li><code>`--pp`</code> - path to individual JSON file corresponding to Phenopacket</li>
 * <li><code>`--output-scores`</code> - path where results in TSV format will be written</li>
 * </ul>
 * </p>
 */
@Component
public class ScorePhenopacketsRunner implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(ScorePhenopacketsRunner.class);

    private static final String DELIMITER = "\t";

    private final GenomeSequenceAccessor genomeSequenceAccessor;

    private final SplicingTranscriptSource splicingTranscriptSource;

    private final SplicingEvaluator splicingEvaluator;

    /**
     * List of paths to phenopacket JSONs which will be analyzed.
     */
    private final List<Path> phenopacketPaths = new ArrayList<>();

    /**
     * Path to TSV file where results will be written.
     */
    private Path outputPath;

    public ScorePhenopacketsRunner(GenomeSequenceAccessor genomeSequenceAccessor, SplicingTranscriptSource splicingTranscriptSource, SplicingEvaluator splicingEvaluator) {
        this.genomeSequenceAccessor = genomeSequenceAccessor;
        this.splicingTranscriptSource = splicingTranscriptSource;
        this.splicingEvaluator = splicingEvaluator;
    }

    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("score-phenopackets")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // complaints raised in the parse cli function
            return;
        }

        LOGGER.info("Analyzing {} phenopackets", phenopacketPaths.size());

        try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
            List<ScoringStrategy> scoringStrategies = Arrays.asList(ScoringStrategy.CANONICAL_DONOR, ScoringStrategy.CRYPTIC_DONOR, ScoringStrategy.CRYPTIC_DONOR_IN_CANONICAL_POSITION,
                    ScoringStrategy.CANONICAL_ACCEPTOR, ScoringStrategy.CRYPTIC_ACCEPTOR, ScoringStrategy.CRYPTIC_ACCEPTOR_IN_CANONICAL_POSITION);

            List<String> headerFields = new ArrayList<>(Arrays.asList("PHENOPACKET", "VARIANT", "TRANSCRIPT", "MAX_SCORE"));
            for (ScoringStrategy strategy : scoringStrategies) {
                headerFields.add(strategy.toString());
            }

            // write header
            writer.write(String.join(DELIMITER, headerFields));
            writer.newLine();

            for (Path phenopacketPath : phenopacketPaths) {
                Phenopacket phenopacket = Utils.readPhenopacket(phenopacketPath);
                String phenopacketName = phenopacketPath.toFile().getName();

                for (Variant variant : phenopacket.getVariantsList()) {
                    if (!variant.getAlleleCase().equals(Variant.AlleleCase.VCF_ALLELE)) {
                        LOGGER.info("Variant allele is not in VCF format: {}\nSkipping..", variant);
                        continue;
                    }
                    // make variant proper for splicing analysis
                    VcfAllele vcfAllele = variant.getVcfAllele();
                    GenomeCoordinates varCoordinates = GenomeCoordinates.newBuilder()
                            .setContig(vcfAllele.getChr())
                            .setBegin(vcfAllele.getPos() - 1)
                            .setEnd(vcfAllele.getPos() + vcfAllele.getRef().length() - 1)
                            .setStrand(true)
                            .build();

                    SplicingVariant splv = SplicingVariant.newBuilder()
                            .setCoordinates(varCoordinates)
                            .setRef(vcfAllele.getRef())
                            .setAlt(vcfAllele.getAlt())
                            .build();
                    LOGGER.info("Evaluating variant {}", splv);

                    // fetch nucleotide sequence neighboring the variant
                    SequenceInterval sequenceInterval = genomeSequenceAccessor.fetchSequence(varCoordinates.getContig(),
                            varCoordinates.getBegin() - 50, varCoordinates.getEnd() + 50,
                            true);

                    // fetch all the transcripts overlapping with variant's position
                    List<SplicingTranscript> transcripts = splicingTranscriptSource.fetchTranscripts(varCoordinates.getContig(), varCoordinates.getBegin(), varCoordinates.getEnd());

                    // evaluate variant against all the transcripts & write out the scores
                    for (SplicingTranscript transcript : transcripts) {
                        SplicingPathogenicityData evaluation = splicingEvaluator.evaluate(splv, transcript, sequenceInterval);
                        StringBuilder builder = new StringBuilder()
                                .append(phenopacketName).append(DELIMITER)
                                .append(String.format("%s:%d %s>%s", splv.getContig(), splv.getPos(), splv.getRef(), splv.getAlt())).append(DELIMITER)
                                .append(transcript.getAccessionId()).append(DELIMITER)
                                .append(evaluation.getMaxScore()); // no delimiter here!

                        // write out all the scores
                        for (ScoringStrategy strategy : scoringStrategies) {
                            builder.append(DELIMITER).append(evaluation.getScoresMap().getOrDefault(strategy, Double.NaN));
                        }

                        writer.write(builder.toString());
                        writer.newLine();
                    }
                }
            }
        }

    }

    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopackets - either all from a directory, or specified by `--pp` option
        noJsonFiles:
        if (args.containsOption("pp-dir")) {
            String ppDirString = args.getOptionValues("pp-dir").get(0);
            Path ppDirPath = Paths.get(ppDirString);
            File[] jsonFiles = ppDirPath.toFile().listFiles(f -> f.getName().endsWith(".json"));
            if (jsonFiles == null) {
                break noJsonFiles;
            }
            Arrays.stream(jsonFiles).map(File::toPath).forEach(phenopacketPaths::add);
        }

        List<String> pps = args.getOptionValues("pp");
        if (pps != null) {
            phenopacketPaths.addAll(pps.stream().map(Paths::get).collect(Collectors.toList()));
        }

        // Output file path - where to write TSV file with scores
        if (!args.containsOption("output-scores")) {
            LOGGER.error("Missing '--output-scores' argument");
            return false;
        }
        outputPath = Paths.get(args.getOptionValues("output-scores").get(0));

        return true;
    }
}
