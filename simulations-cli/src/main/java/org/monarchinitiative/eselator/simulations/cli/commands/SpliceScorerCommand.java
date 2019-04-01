package org.monarchinitiative.eselator.simulations.cli.commands;

import de.charite.compbio.jannovar.data.JannovarData;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.monarchinitiative.eselator.simulations.cli.CustomVariantEvaluation;
import org.monarchinitiative.eselator.simulations.cli.Utils;
import org.monarchinitiative.exomiser.core.genome.GenomeAnalysisService;
import org.monarchinitiative.exomiser.core.genome.VariantAnnotator;
import org.monarchinitiative.exomiser.core.genome.dao.splicing.SpliceScorer;
import org.monarchinitiative.exomiser.core.genome.dao.splicing.SplicingContext;
import org.monarchinitiative.exomiser.core.genome.dao.splicing.SplicingDao;
import org.monarchinitiative.exomiser.core.model.VariantAnnotation;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

@Component
public class SpliceScorerCommand implements ApplicationRunner {


    private static final Logger LOGGER = LoggerFactory.getLogger(SpliceScorerCommand.class);

    private final List<String> phenopacketPaths;

    private final GenomeAnalysisService genomeAnalysisService;

    private final VariantAnnotator variantAnnotator;

    private final JannovarData jannovarData;

    private final SplicingDao splicingDao;

    private SpliceScorer scorer;

    private Path outputPath;

    private String strategy;

    public SpliceScorerCommand(GenomeAnalysisService genomeAnalysisService, VariantAnnotator variantAnnotator, JannovarData jannovarData, SplicingDao splicingDao) {
        this.genomeAnalysisService = genomeAnalysisService;
        this.variantAnnotator = variantAnnotator;
        this.jannovarData = jannovarData;
        this.splicingDao = splicingDao;
        this.phenopacketPaths = new ArrayList<>();
    }

    /**
     * @return function for splitting input string such as <code>key=value</code> into an <code>['key', 'value']</code>
     * array. The function is guaranteed to return array with at least 2 fields. For string <code>key=</code>, the array
     * will look like <code>['key', '']</code>
     */
    private static Function<String, String[]> decodeInfoField() {
        return s -> {
            if (s.contains("=")) {
                return s.split("=", 2);
            } else {
                return new String[]{s, ""};
            }
        };
    }

    private static Object[] makeRecordRepresentation(CustomVariantEvaluation cve, Map<String, String> fields) {
        VariantEvaluation ve = cve.getVariantEvaluation();
        return new Object[]{
                // PP_ID
                cve.getPhenopacket().getId(),
                // VARIANT
                String.format("%s:%d%s>%s", ve.getChromosomeName(), ve.getPosition(), ve.getRef(), ve.getAlt()),
                // VCLASS
                fields.get("VCLASS"),
                // PATHOMECHANISM
                fields.get("PATHOMECHANISM"),
                // CONSEQUENCE
                fields.get("CONSEQUENCE"),
                // DONOR
                cve.getDonorScore(),
                // ACCEPTOR
                cve.getAcceptorScore(),
                // EXON
                cve.getExonScore(),
                // INTRON
                cve.getIntronScore()};
    }

    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("splice-scorer")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }

        // ----------------- SELECT SCORER ------------------------------------
        if (!splicingDao.getScorersMap().containsKey(strategy)) {
            LOGGER.warn("Unknown strategy '{}'. Known strategies: '{}'", strategy, splicingDao.getScorersMap().keySet());
            return;
        } else {
            scorer = splicingDao.getScorersMap().get(strategy);
        }
        LOGGER.info("Using scoring strategy '{}'", strategy);

        // ----------------- READ PHENOPACKETS --------------------------------
        List<Phenopacket> ppList = new ArrayList<>();
        for (String ppath : phenopacketPaths) {
            Path p = Paths.get(ppath);
            try {
                final Phenopacket pp = Utils.readPhenopacket(p);
                ppList.add(pp);
            } catch (IOException ioe) {
                LOGGER.warn("Unable to read phenopacket from '{}'. ", ppath, ioe);
            }
        }
        LOGGER.info("Read {} phenopackets", ppList.size());

        // ----------------- SCORE VARIANTS -----------------------------------
        LOGGER.info("Scoring variants");
        // map to custom variant evaluations first
        List<CustomVariantEvaluation> cveList = new ArrayList<>();
        for (Phenopacket pp : ppList) {
            for (Variant pv : pp.getVariantsList()) {
                switch (pv.getAlleleCase()) {
                    case HGVS_ALLELE:
                    case ISCN_ALLELE:
                    case SPDI_ALLELE:
                    case ALLELE_NOT_SET:
                    case MOUSE_ALLELE:
                        LOGGER.warn("Allele format '{}' not supported at the moment (pp: '{}', var: '{}'", pv.getAlleleCase(), pp.getId(), pv);
                        continue;
                    case VCF_ALLELE:
                        VcfAllele va = pv.getVcfAllele();
                        Integer chr = jannovarData.getRefDict().getContigNameToID().get(va.getChr());
                        VariantAnnotation annotation = variantAnnotator.annotate(va.getChr(), va.getPos(), va.getRef(), va.getAlt());
                        VariantEvaluation ve = VariantEvaluation.builder(chr, va.getPos(), va.getRef(), va.getAlt())
                                .annotations(annotation.getTranscriptAnnotations())
                                .build();

                        SplicingContext sctx = splicingDao.buildSplicingContext(ve);
                        CustomVariantEvaluation cv = CustomVariantEvaluation.builder()
                                .setPhenopacket(pp)
                                .setVariantEvaluation(ve)
                                .setVcfAlleleInfoField(va.getInfo())
                                .setDonorScore(scorer.pathogenicityForDonor(sctx))
                                .setAcceptorScore(scorer.pathogenicityForAcceptor(sctx))
                                .setIntronScore(scorer.pathogenicityForIntron(sctx))
                                .setExonScore(scorer.pathogenicityForExon(sctx))
                                .build();
                        cveList.add(cv);
                }
            }
        }

        // ----------------- WRITE CUSTOM VARIANT EVALUATIONS TO FILE ---------
        try (CSVPrinter printer = new CSVPrinter(Files.newBufferedWriter(outputPath),
                CSVFormat.TDF
                        .withHeader("PP_ID", "VARIANT", "VCLASS", "PATHOMECHANISM", "CONSEQUENCE", "DONOR", "ACCEPTOR", "EXON", "INTRON"))) {

            for (CustomVariantEvaluation cve : cveList) {
                // Decode INFO string:
                // token[0] = id (e.g. VCLASS), token[1] = value (e.g. coding)
                Map<String, String> variantInfoFields = Arrays.stream(cve.getVcfAlleleInfoField().split(";"))
                        .map(decodeInfoField())
                        .collect(Collectors.toMap(s -> s[0], s -> s[1]));

                // write the variant
                printer.printRecord(makeRecordRepresentation(cve, variantInfoFields));
            }
        }
        LOGGER.info("Wrote {} variants", cveList.size());
        LOGGER.info("Done!");
    }

    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopackets
        if (!args.containsOption("pp")) {
            LOGGER.warn("At least one '--pp' argument for Phenopacket path must be present");
            return false;
        }
        phenopacketPaths.addAll(args.getOptionValues("pp"));

        // Output file path - results
        if (!args.containsOption("output")) {
            LOGGER.warn("Missing '--output' argument");
            return false;
        }
        outputPath = Paths.get(args.getOptionValues("output").get(0));

        if (!args.containsOption("strategy")) {
            LOGGER.warn("Missing '--strategy' argument");
        }
        strategy = args.getOptionValues("strategy").get(0);

        return true;
    }
}
