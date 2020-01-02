package org.monarchinitiative.threes.benchmarks.cmd;

import de.charite.compbio.jannovar.data.ReferenceDictionary;
import de.charite.compbio.jannovar.reference.*;
import net.sourceforge.argparse4j.inf.MutuallyExclusiveGroup;
import net.sourceforge.argparse4j.inf.Namespace;
import net.sourceforge.argparse4j.inf.Subparser;
import net.sourceforge.argparse4j.inf.Subparsers;
import org.monarchinitiative.threes.benchmarks.Utils;
import org.monarchinitiative.threes.core.data.SplicingTranscriptSource;
import org.monarchinitiative.threes.core.model.SplicingTranscript;
import org.monarchinitiative.threes.core.scoring.SplicingAnnotator;
import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;
import org.monarchinitiative.threes.core.scoring.sparse.ScoringStrategy;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;
import xyz.ielis.hyperutil.reference.fasta.GenomeSequenceAccessor;
import xyz.ielis.hyperutil.reference.fasta.SequenceInterval;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
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
public class ScorePhenopacketsCommand extends Command {

    private static final Logger LOGGER = LoggerFactory.getLogger(ScorePhenopacketsCommand.class);

    private static final String DELIMITER = "\t";

    private final GenomeSequenceAccessor genomeSequenceAccessor;

    private final SplicingTranscriptSource splicingTranscriptSource;

    private final SplicingAnnotator splicingAnnotator;

    /**
     * List of paths to phenopacket JSONs which will be analyzed.
     */
    private final List<Path> phenopacketPaths = new ArrayList<>();

    public ScorePhenopacketsCommand(GenomeSequenceAccessor genomeSequenceAccessor,
                                    SplicingTranscriptSource splicingTranscriptSource,
                                    SplicingAnnotator splicingAnnotator) {
        this.genomeSequenceAccessor = genomeSequenceAccessor;
        this.splicingTranscriptSource = splicingTranscriptSource;
        this.splicingAnnotator = splicingAnnotator;
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

    public static void setupSubparsers(Subparsers subparsers) {
        Subparser ic = subparsers.addParser("score-phenopackets")
                .setDefault("cmd", "score-phenopackets")
                .help("score phenopackets");

        MutuallyExclusiveGroup ppGroup = ic.addMutuallyExclusiveGroup();
        ppGroup.addArgument("--pp-dir")
                .help("path to ClinVar VCF file");

        ppGroup.addArgument("--pp")
                .nargs("*")
                .help("path");

        ic.addArgument("output-scores")
                .help("where to write the output file");
    }

    @Override
    public void run(Namespace namespace) throws CommandException {
        // 0 - parse CLI args
        // Collect paths from directory, or individual paths specified by `--pp` option
        final String ppDir = namespace.getString("pp_dir");
        if (ppDir != null) {
            final Path ppDirPath = Paths.get(ppDir);
            LOGGER.info("Reading phenopackets from '{}'", ppDirPath);

            File[] jsonFiles = ppDirPath.toFile().listFiles(f -> f.getName().endsWith(".json"));
            if (jsonFiles != null) {
                Arrays.stream(jsonFiles).map(File::toPath).forEach(phenopacketPaths::add);
            }
        }

        final List<String> pps = namespace.getList("pp");
        if (pps != null) {
            phenopacketPaths.addAll(pps.stream().map(Paths::get).collect(Collectors.toList()));
        }

        /*
         * Path to TSV file where results will be written.
         */
        Path outputPath = Paths.get(namespace.getString("output_scores"));


        // 1 - run the app

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

            List<String> headerFields = new ArrayList<>(Arrays.asList("CASE",
                    "VARIANT", "TRANSCRIPT", "VCLASS", "PATHOMECHANISM", "CONSEQUENCE",
                    "MAX_SCORE"));
            for (ScoringStrategy strategy : scoringStrategies) {
                headerFields.add(strategy.toString());
            }

            // write header
            writer.write(String.join(DELIMITER, headerFields));
            writer.newLine();

            // analyze phenopackets
            final ReferenceDictionary rd = genomeSequenceAccessor.getReferenceDictionary();
            LOGGER.info("Analyzing {} phenopackets", phenopacketPaths.size());
            for (Path phenopacketPath : phenopacketPaths) {
                Phenopacket phenopacket = Utils.readPhenopacket(phenopacketPath);

                for (Variant variant : phenopacket.getVariantsList()) {
                    if (!variant.getAlleleCase().equals(Variant.AlleleCase.VCF_ALLELE)) {
                        LOGGER.info("Variant allele is not in VCF format: {}\nSkipping..", variant);
                        continue;
                    }


                    // make variant proper for splicing analysis
                    VcfAllele vcfAllele = variant.getVcfAllele();
                    if (!rd.getContigNameToID().containsKey(vcfAllele.getChr())) {
                        LOGGER.info("Unknown chromosome of variant {}:{}{}>{}", vcfAllele.getChr(), vcfAllele.getPos(), vcfAllele.getRef(), vcfAllele.getAlt());
                        continue;
                    }

                    final GenomeVariant gv = new GenomeVariant(new GenomePosition(rd, Strand.FWD, rd.getContigNameToID().get(vcfAllele.getChr()), vcfAllele.getPos(), PositionType.ONE_BASED),
                            vcfAllele.getRef(), vcfAllele.getAlt());

                    LOGGER.debug("Evaluating variant `{}`", gv);


                    // fetch all the transcripts overlapping with variant's position
                    // retain the curated transcripts (accession id starts with 'NM_')
                    List<SplicingTranscript> curatedRefSeqTranscripts = splicingTranscriptSource.fetchTranscripts(gv.getChrName(), gv.getGenomeInterval().getBeginPos(), gv.getGenomeInterval().getEndPos(), rd).stream()
                            .filter(tx -> tx.getAccessionId().startsWith("NM_"))
                            .collect(Collectors.toList());

                    if (curatedRefSeqTranscripts.isEmpty()) {
                        LOGGER.warn("No curated RefSeq transcript overlaps with variant `{}`", gv);
                        continue;
                    }


                    SplicingTranscript transcript = curatedRefSeqTranscripts.stream()
                            .min(Utils.transcriptPriorityComparator())
                            .get();


                    // fetch nucleotide sequence neighboring the variant
                    final GenomeInterval txQueryRegion = transcript.getTxRegionCoordinates().withMorePadding(50, 50);
                    Optional<SequenceInterval> sequenceIntervalOpt = genomeSequenceAccessor.fetchSequence(txQueryRegion);
                    if (!sequenceIntervalOpt.isPresent()) {
                        LOGGER.warn("Unable to fetch sequence `{}` for transcript `{}-{}`", txQueryRegion, transcript.getAccessionId(), transcript.getTxRegionCoordinates());
                        continue;
                    }

                    // get VCLASS, PATHOMECHANISM, CONSEQUENCE
                    Map<String, String> infos = getInfoFromVcfAllele(vcfAllele.getInfo());

                    // --- EVALUATE VARIANT AGAINST THE TRANSCRIPT & WRITE OUT THE SCORES ---
                    // evaluate
                    SplicingPathogenicityData evaluation = splicingAnnotator.evaluate(gv, transcript, sequenceIntervalOpt.get());

                    // write out
                    StringBuilder builder = new StringBuilder()
                            .append(phenopacket.getId()).append(DELIMITER)
                            // VARIANT
                            .append(String.format("%s:%d %s>%s", gv.getChrName(), gv.getPos() + 1, gv.getRef(), gv.getAlt())).append(DELIMITER)
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
                }
            }
        } catch (IOException e) {
            LOGGER.error("Error: ", e);
            throw new CommandException(e);
        }

        LOGGER.info("くまくま━━━━━━ヽ（ ・(ｪ)・ ）ノ━━━━━━ !!!");
        LOGGER.info("                 Done!               ");
    }
}