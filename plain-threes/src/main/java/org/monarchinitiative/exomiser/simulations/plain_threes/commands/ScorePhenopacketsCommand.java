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
import java.util.*;
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
 * <b>!! IMPORTANT !!</b> - this code does not work with other than RefSeq splicing transcript source.
 */
@Component
public class ScorePhenopacketsCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(ScorePhenopacketsCommand.class);

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

    public ScorePhenopacketsCommand(GenomeSequenceAccessor genomeSequenceAccessor,
                                    SplicingTranscriptSource splicingTranscriptSource,
                                    SplicingEvaluator splicingEvaluator) {
        this.genomeSequenceAccessor = genomeSequenceAccessor;
        this.splicingTranscriptSource = splicingTranscriptSource;
        this.splicingEvaluator = splicingEvaluator;
    }

    private static Map<String, String> getInfoFromVcfAllele(String info) {
        Map<String, String> map = new HashMap<>();
        // info looks like 'VCLASS=splicing;PATHOMECHANISM=splicing|5ss|disrupted;CONSEQUENCE=Exon skipping'
        String[] field = info.split(";");
        for (String tokens : field) {
            // tokens look like 'VCLASS=splicing', ...

            String[] token = tokens.split("=");
            if (token.length == 1) {
                map.put(token[0], "None");
            } else if (token.length == 2) {
                map.put(token[0], token[1]);
            }
        }

        return map;
    }

    /**
     * @return {@link Comparator} where the highest priority transcript is the one with the smallest integer value of the
     * central part of the accession ID. The central part of the <em>NM_004004.2</em> is <em>004004</em>.
     */
    private static Comparator<? super SplicingTranscript> transcriptPriorityComparator() {
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


        try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
            // prepare header
            List<ScoringStrategy> scoringStrategies = Arrays.asList(
                    // DONOR
                    ScoringStrategy.CANONICAL_DONOR, ScoringStrategy.CRYPTIC_DONOR, ScoringStrategy.CRYPTIC_DONOR_IN_CANONICAL_POSITION,
                    // ACCEPTOR
                    ScoringStrategy.CANONICAL_ACCEPTOR, ScoringStrategy.CRYPTIC_ACCEPTOR, ScoringStrategy.CRYPTIC_ACCEPTOR_IN_CANONICAL_POSITION,
                    // ESE/ESS
                    ScoringStrategy.SMS
            );

            List<String> headerFields = new ArrayList<>(Arrays.asList("PHENOPACKET",
                    "VARIANT", "TRANSCRIPT", "VCLASS", "PATHOMECHANISM", "CONSEQUENCE",
                    "MAX_SCORE"));
            for (ScoringStrategy strategy : scoringStrategies) {
                headerFields.add(strategy.toString());
            }

            // write header
            writer.write(String.join(DELIMITER, headerFields));
            writer.newLine();

            // analyze phenopackets
            LOGGER.info("Analyzing {} phenopackets", phenopacketPaths.size());
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
                    LOGGER.debug("Evaluating variant {}", splv);


                    // fetch all the transcripts overlapping with variant's position
                    // retain the curated transcripts (accession id starts with 'NM_')
                    List<SplicingTranscript> curatedTranscripts = splicingTranscriptSource.fetchTranscripts(varCoordinates.getContig(), varCoordinates.getBegin(), varCoordinates.getEnd()).stream()
                            .filter(tx -> tx.getAccessionId().startsWith("NM_"))
                            .collect(Collectors.toList());

                    if (curatedTranscripts.isEmpty()) {
                        LOGGER.warn("No curated transcripts overlap with variant {}", splv);
                        break;
                    }


                    SplicingTranscript transcript = curatedTranscripts.stream()
                            .min(transcriptPriorityComparator())
                            .get();


                    // fetch nucleotide sequence neighboring the variant
                    SequenceInterval sequenceInterval = genomeSequenceAccessor.fetchSequence(transcript.getContig(),
                            transcript.getTxBegin() - 50, transcript.getTxEnd() + 50,
                            transcript.getStrand());

                    // get VCLASS, PATHOMECHANISM, CONSEQUENCE
                    Map<String, String> infos = getInfoFromVcfAllele(vcfAllele.getInfo());

                    // evaluate variant against all the transcripts & write out the scores
//                    for (SplicingTranscript transcript : curatedTranscripts) {

                    // evaluate
                    SplicingPathogenicityData evaluation = splicingEvaluator.evaluate(splv, transcript, sequenceInterval);

                    // write out
                    StringBuilder builder = new StringBuilder()
                            .append(phenopacketName).append(DELIMITER)
                            // VARIANT
                            .append(String.format("%s:%d %s>%s", splv.getContig(), splv.getPos(), splv.getRef(), splv.getAlt())).append(DELIMITER)
                            // TRANSCRIPT
                            .append(transcript.getAccessionId()).append(DELIMITER)
                            // VCLASS
                            .append(infos.getOrDefault("VCLASS", "None")).append(DELIMITER)
                            // PATHOMECHANISM
                            .append(infos.getOrDefault("PATHOMECHANISM", "None")).append(DELIMITER)
                            // CONSEQUENCE
                            .append(infos.getOrDefault("CONSEQUENCE", "None")).append(DELIMITER)
                            // MAX SCORE
                            .append(evaluation.getMaxScore()); // no delimiter here!

                    // write out all the scores
                    for (ScoringStrategy strategy : scoringStrategies) {
                        builder.append(DELIMITER).append(evaluation.getScoresMap().getOrDefault(strategy, Double.NaN));
                    }

                    writer.write(builder.toString());
                    writer.newLine();
//                    }
                }
            }
        }

        LOGGER.info("くまくま━━━━━━ヽ（ ・(ｪ)・ ）ノ━━━━━━ !!!");
        LOGGER.info("                 Done!               ");
    }

    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopackets - either all from a directory, or specified by `--pp` option
        noJsonFiles:
        if (args.containsOption("pp-dir")) {
            String ppDirString = args.getOptionValues("pp-dir").get(0);
            Path ppDirPath = Paths.get(ppDirString);
            LOGGER.info("Reading phenopackets from '{}'", ppDirPath);
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
