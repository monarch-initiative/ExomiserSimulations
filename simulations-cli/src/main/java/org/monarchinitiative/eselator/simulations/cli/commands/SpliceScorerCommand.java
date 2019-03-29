package org.monarchinitiative.eselator.simulations.cli.commands;

import de.charite.compbio.jannovar.data.JannovarData;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVPrinter;
import org.monarchinitiative.eselator.simulations.cli.CustomVariantEvaluation;
import org.monarchinitiative.eselator.simulations.cli.Utils;
import org.monarchinitiative.exomiser.core.genome.GenomeAnalysisService;
import org.monarchinitiative.exomiser.core.genome.dao.splicing.SpliceScorer;
import org.monarchinitiative.exomiser.core.genome.dao.splicing.SplicingContext;
import org.monarchinitiative.exomiser.core.genome.dao.splicing.SplicingDao;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.Variant;
import org.phenopackets.schema.v1.core.VcfAllele;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

@Component
public class SpliceScorerCommand implements ApplicationRunner {


    private static final Logger LOGGER = LoggerFactory.getLogger(SpliceScorerCommand.class);

    private final List<String> phenopacketPaths;

    private final GenomeAnalysisService genomeAnalysisService;

    private final JannovarData jannovarData;

    private final SplicingDao splicingDao;

    private SpliceScorer scorer;

    private Path outputPath;

    private String strategy;

    public SpliceScorerCommand(GenomeAnalysisService genomeAnalysisService, JannovarData jannovarData, SplicingDao splicingDao) {
        this.genomeAnalysisService = genomeAnalysisService;
        this.jannovarData = jannovarData;
        this.splicingDao = splicingDao;
        this.phenopacketPaths = new ArrayList<>();
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

        // ----------------- MAP PHENOPACKETS TO CUSTOM VARIANT EVALUATIONS ---
        LOGGER.info("Scoring variants");
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
                        VariantEvaluation ve = VariantEvaluation.builder(chr, va.getPos(), va.getRef(), va.getAlt()).build();
                        SplicingContext sctx = splicingDao.buildSplicingContext(ve);
                        CustomVariantEvaluation cv = CustomVariantEvaluation.builder()
                                .setPhenopacket(pp)
                                .setVariantEvaluation(ve)
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
                        .withHeader("PP_ID", "VARIANT", "DONOR", "ACCEPTOR", "EXON", "INTRON"))) {

            for (CustomVariantEvaluation cve : cveList) {
                printer.printRecord(
                        cve.getPhenopacket().getId(), // PP_ID
                        cve.getVariantEvaluation(), // VARIANT
                        cve.getDonorScore(), // DONOR
                        cve.getAcceptorScore(), // ACCEPTOR
                        cve.getExonScore(), // EXON
                        cve.getIntronScore() // INTRON
                        );
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
