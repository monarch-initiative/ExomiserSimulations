package org.monarchinitiative.exomiser_simulations.cli.commands;

import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.model.GeneScore;
import org.monarchinitiative.exomiser_simulations.cli.Utils;
import org.monarchinitiative.exomiser_simulations.cli.simulators.SingleVcfSimulator;
import org.monarchinitiative.exomiser_simulations.cli.simulators.VcfSimulator;
import org.monarchinitiative.exomiser.core.Exomiser;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisMode;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.phenopackets.schema.v1.Phenopacket;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static java.util.Map.Entry.comparingByKey;

@Component
public class LiricalCommand implements ApplicationRunner {
    private static final Logger LOGGER = LoggerFactory.getLogger(LiricalCommand.class);

    /**
     * HARDCODED FOR NOW
     */
    private static final float MAX_FREQ = 0.001f;

    private final Exomiser exomiser;

    private Path templateVcfPath;

    private Path phenopacketDirectoryPath;

    private List<String> resultlist=new ArrayList<>();
    /** key -- a rank, e.g., 3; value -- number of diseases with this rank inthe simulation. */
    private Map<Integer,Integer> rankCountMap;



    public LiricalCommand(Exomiser exomiser) {
        this.exomiser = exomiser;
        this.rankCountMap = new HashMap<>();
    }

    /**
     * Only present (non-negated) PhenotypicFeatures are reported
     *
     * @param pp {@link Phenopacket} describing the proband
     * @return list of HPO id strings representing subject's phenotype
     */
    private static List<String> getPresentPhenotypesAsHpoStrings(Phenopacket pp) {
        return pp.getPhenotypicFeaturesList().stream()
                .filter(p -> !p.getNegated())
                .map(p -> p.getType().getId())
                .collect(Collectors.toList());
    }

    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (!args.containsOption("lirical")) {
            // not running this command
            return;
        }

        if (!parseCliArgs(args)) {
            // unable to parse command line, complaints raised in the function
            return;
        }

        // -----------------------    FOR EACH PHENOPACKET    --------------------------------------
        File[] fileArray = phenopacketDirectoryPath.toFile().listFiles();
        if (fileArray == null) {
            LOGGER.warn("Phenopacket file array is null. This should not happen");
            return;
        }
        for (File phenopacketFilePath : fileArray) {
            // -----------------------    READ PHENOPACKET    --------------------------------------
            LOGGER.info("Reading phenopacket from '{}'", phenopacketFilePath);
            Phenopacket pp;
            pp = Utils.readPhenopacket(phenopacketFilePath.toPath());
            if (pp.getSubject().getId().isEmpty()) {
                LOGGER.warn("Phenopacket subject's ID must not be empty. Unable to continue");
                System.exit(1);
            }

            if (pp.getGenesCount() != 1) {
                LOGGER.error("[ERROR] Phenopackets used for simulation MUST have exactly one gene");
                System.exit(1);
            }
            String entrezString = pp.getGenes(0).getId();
            if (entrezString.contains("NCBIGene:")) {
                entrezString = entrezString.substring(9);
            }
            int entrezId = Integer.parseInt(entrezString);

            String symbol = pp.getGenes(0).getSymbol();
            if (pp.getDiseasesCount()!=1) {
                LOGGER.error("[ERROR] Phenopackets used for simulation MUST have exactly one disease");
                System.exit(1);
            }
            String diseaseId = pp.getDiseases(0).getTerm().getId();
            String diseaseLabel = pp.getDiseases(0).getTerm().getLabel();

            // -----------------------    CREATE THE SIMULATED VCF FILE    -------------------------
            LOGGER.info("Creating simulated VCF file");
            VcfSimulator simulator = new SingleVcfSimulator(templateVcfPath);
            Path vcfPath = simulator.simulateVcfWithPhenopacket(pp);


            // -----------------------    FORGE EXOMISER ANALYSIS    -------------------------------
            Set<FrequencySource> frequencySources = new HashSet<>(FrequencySource.ALL_EXTERNAL_FREQ_SOURCES);
            LOGGER.info("Creating analysis");
            Analysis analysis = exomiser.getAnalysisBuilder()
                    .analysisMode(AnalysisMode.PASS_ONLY)
                    .vcfPath(vcfPath)
                    .probandSampleName(pp.getSubject().getId())
                    .genomeAssembly(GenomeAssembly.HG19)
                    .frequencySources(frequencySources)
                    .addFrequencyFilter(MAX_FREQ)
                    .hpoIds(getPresentPhenotypesAsHpoStrings(pp))
                    .addOmimPrioritiser()
                    .addHiPhivePrioritiser()
                    .build();


            // -----------------------    RUN THE ANALYSIS AND WRITE THE RESULTS    ----------------
            LOGGER.info("Running the analysis");
            final AnalysisResults results = exomiser.run(analysis);
            int rank = getRankOfGene(results,entrezId);
            // this is the rank of the target gene
            String res = String.format("[INFO] Rank of gene %s [%s;%s] was %d", symbol,diseaseLabel,diseaseId,rank);
            System.out.println(res);
            resultlist.add(res);
            rankCountMap.putIfAbsent(rank,0);
            int c = 1 + rankCountMap.get(rank);
            rankCountMap.put(rank,c);
            System.out.println(results);
        }
        printOutSimulationResults();
    }



    private int getRankOfGene(AnalysisResults ares,int targetEntrezGeneId) {
        List<Gene> genescores = ares.getGenes();
        int N = genescores.size();
        int rank = 0;
        for (Gene gene :genescores) {
            int entrez = gene.getEntrezGeneID();
            String symb = gene.getGeneSymbol();
            rank++;
            if (entrez == targetEntrezGeneId) {
                // this is the rank of the target gene
               return rank;
            }
        }
        return N+1; // tied for the rank above the worst rank of
        // all genes that were actually ranked.
    }


    private void printOutSimulationResults() throws  IOException{
        BufferedWriter writer = new BufferedWriter(new FileWriter("exomiser_sim_results.txt"));

        // sort the map by values first
        Map<Integer, Integer> sorted = this.rankCountMap
                .entrySet()
                .stream()
                //.sorted(comparingByValue())
                .sorted(comparingByKey())
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2,
                        LinkedHashMap::new));
        int total = this.rankCountMap.values().stream().mapToInt(i->i).sum();
        for (Map.Entry<Integer, Integer> e : sorted.entrySet()) {
            System.out.println(String.format("%s: %d (%.1f%%)", e.getKey(), e.getValue(),(100.0*e.getValue()/total)));
            writer.write(String.format("%s: %d (%.1f%%)", e.getKey(), e.getValue(),(100.0*e.getValue()/total)) + "\n");
        }
        writer.close();
    }






    /**
     * Parse the command line arguments and return false if anything is missing. Errors are also logged.
     */
    private boolean parseCliArgs(ApplicationArguments args) {
        // Phenopacket
        if (!args.containsOption("pp")) {
            LOGGER.warn("Missing '--pp' argument for Phenopacket directory path");
            return false;
        }
        phenopacketDirectoryPath = Paths.get(args.getOptionValues("pp").get(0));
        if (!phenopacketDirectoryPath.toFile().isDirectory()) {
            LOGGER.warn("Path '{}' does not point to existing directory", phenopacketDirectoryPath);
            return false;
        }

        // Template VCF file
        if (!args.containsOption("vcf")) {
            LOGGER.warn("Missing '--vcf' argument for template VCF file path");
            return false;
        }
        templateVcfPath = Paths.get(args.getOptionValues("vcf").get(0));

        return true;
    }
}
