package org.monarchinitiative.exomiser.simulations.cli.commands;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import de.charite.compbio.jannovar.mendel.SubModeOfInheritance;
import org.monarchinitiative.exomiser.core.Exomiser;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisMode;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.util.InheritanceModeOptions;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
import org.monarchinitiative.exomiser.simulations.cli.Utils;
import org.monarchinitiative.exomiser.simulations.cli.simulators.SingleVcfSimulator;
import org.monarchinitiative.exomiser.simulations.cli.simulators.VcfSimulator;
import org.phenopackets.schema.v1.Phenopacket;
import org.phenopackets.schema.v1.core.PhenotypicFeature;
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

@Component
public class LiricalCommand implements ApplicationRunner {

    private static final Logger LOGGER = LoggerFactory.getLogger(LiricalCommand.class);

    /**
     * HARDCODED FOR NOW
     */
    // Exomiser uses percentage values, so this is maybe too low?
    private static final float MAX_FREQ = 0.001f;

    private static final InheritanceModeOptions INHERITANCE_MODE_OPTIONS;
    static {
        Map<SubModeOfInheritance, Float> inheritanceModeFrequencyCutoffs = new EnumMap<>(SubModeOfInheritance.class);
        // all frequencies are in percentage values
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.AUTOSOMAL_DOMINANT, 0.1f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.AUTOSOMAL_RECESSIVE_COMP_HET, 2.0f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.AUTOSOMAL_RECESSIVE_HOM_ALT, 0.1f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.X_DOMINANT, 0.1f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.X_RECESSIVE_COMP_HET, 2.0f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.X_RECESSIVE_HOM_ALT, 0.1f);
        inheritanceModeFrequencyCutoffs.put(SubModeOfInheritance.MITOCHONDRIAL, 0.2f);
        INHERITANCE_MODE_OPTIONS = InheritanceModeOptions.of(inheritanceModeFrequencyCutoffs);
    }

    private static final EnumSet<FrequencySource> FREQUENCY_SOURCES = EnumSet.of(
            FrequencySource.ESP_AFRICAN_AMERICAN, FrequencySource.ESP_ALL, FrequencySource.ESP_EUROPEAN_AMERICAN,
            FrequencySource.THOUSAND_GENOMES,
            FrequencySource.EXAC_AFRICAN_INC_AFRICAN_AMERICAN, FrequencySource.EXAC_AMERICAN, FrequencySource.EXAC_EAST_ASIAN, FrequencySource.EXAC_FINNISH,
            FrequencySource.EXAC_NON_FINNISH_EUROPEAN, FrequencySource.EXAC_SOUTH_ASIAN, FrequencySource.EXAC_OTHER,
            FrequencySource.UK10K, FrequencySource.TOPMED,
            FrequencySource.GNOMAD_E_AFR, FrequencySource.GNOMAD_E_AMR,
            // including the GNOMAD_E_ASJ and GNOMAD_G_ASJ reduces performance on 100K genomes so we disable it by default
            //FrequencySource.GNOMAD_E_ASJ,
            FrequencySource.GNOMAD_E_EAS, FrequencySource.GNOMAD_E_FIN,
            FrequencySource.GNOMAD_E_NFE, FrequencySource.GNOMAD_E_OTH, FrequencySource.GNOMAD_E_SAS,
            FrequencySource.GNOMAD_G_AFR, FrequencySource.GNOMAD_G_AMR,
            //FrequencySource.GNOMAD_G_ASJ,
            FrequencySource.GNOMAD_G_EAS, FrequencySource.GNOMAD_G_FIN,
            FrequencySource.GNOMAD_G_NFE, FrequencySource.GNOMAD_G_OTH, FrequencySource.GNOMAD_G_SAS);

    private static final EnumSet<VariantEffect> NON_EXONIC_EFFECTS = EnumSet.of(
            VariantEffect.FIVE_PRIME_UTR_EXON_VARIANT,
            VariantEffect.FIVE_PRIME_UTR_INTRON_VARIANT,
            VariantEffect.THREE_PRIME_UTR_EXON_VARIANT,
            VariantEffect.THREE_PRIME_UTR_INTRON_VARIANT,
            VariantEffect.NON_CODING_TRANSCRIPT_EXON_VARIANT,
            VariantEffect.NON_CODING_TRANSCRIPT_INTRON_VARIANT,
            VariantEffect.CODING_TRANSCRIPT_INTRON_VARIANT,
            VariantEffect.UPSTREAM_GENE_VARIANT,
            VariantEffect.DOWNSTREAM_GENE_VARIANT,
            VariantEffect.INTERGENIC_VARIANT,
            VariantEffect.REGULATORY_REGION_VARIANT
    );

    private final Exomiser exomiser;

    private Path templateVcfPath;

    private Path phenopacketDirectoryPath;

    private List<String> resultlist=new ArrayList<>();
    /** key -- a rank, e.g., 3; value -- number of diseases with this rank inthe simulation. */
    private Map<Integer,Integer> rankCountMap = new HashMap<>();



    public LiricalCommand(Exomiser exomiser) {
        this.exomiser = exomiser;
    }

    /**
     * Only present (non-negated) {@link PhenotypicFeature}s are reported
     *
     * @param pp {@link Phenopacket} describing the proband
     * @return list of HPO id strings representing subject's phenotype
     */
    static List<String> getPresentPhenotypesAsHpoStrings(Phenopacket pp) {
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
        List<File> phenopackets = Arrays.asList(fileArray);
        for (File phenopacketFilePath : phenopackets) {
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
            if (entrezString.indexOf("ENTREZ:") >= 0) {
                entrezString = entrezString.substring(7);
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
            LOGGER.info("Creating analysis");
            // This analysis is optimised for *exomes* if running a genome sample we need a different configuration.
            Analysis analysis = exomiser.getAnalysisBuilder()
                    .genomeAssembly(GenomeAssembly.HG19)
                    .vcfPath(vcfPath)
                    .probandSampleName(pp.getSubject().getId())
                    .hpoIds(getPresentPhenotypesAsHpoStrings(pp))
                    .analysisMode(AnalysisMode.PASS_ONLY)
                    .inheritanceModes(INHERITANCE_MODE_OPTIONS)
                    .frequencySources(FREQUENCY_SOURCES)
                    .pathogenicitySources(EnumSet.of(PathogenicitySource.REVEL, PathogenicitySource.MVP))
                    // adds an exome mask
                    .addVariantEffectFilter(NON_EXONIC_EFFECTS)
                    .addFailedVariantFilter()
                    // frequency filter max will be automatically derived from the inheritance mode options
                    .addFrequencyFilter()
                    .addPathogenicityFilter(true)
                    .addInheritanceFilter()
                    .addOmimPrioritiser()
                    .addHiPhivePrioritiser()
                    .build();


            // -----------------------    RUN THE ANALYSIS AND WRITE THE RESULTS    ----------------
            LOGGER.info("Running the analysis");
            final AnalysisResults results = exomiser.run(analysis);
            List<Gene> genescores = results.getGenes();
            int rank = 0;
            for (Gene gene : genescores) {
                int entrez = gene.getEntrezGeneID();
                String symb = gene.getGeneSymbol();
                rank++;
                if (entrez == entrezId) {
                    // this is the rank of the target gene
                    String res = String.format("[INFO] Rank of gene %s [%s;%s] was %d", symb,diseaseLabel,diseaseId,rank);
                    System.out.println(res);
                    resultlist.add(res);
                    rankCountMap.putIfAbsent(rank,0);
                    int c = 1 + rankCountMap.get(rank);
                    rankCountMap.put(rank,c);
                }
            }

            System.out.println(results);
        }
        printOutSimulationResults();
    }


    private void printOutSimulationResults() {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("exomiser_sim_results.txt"));
            int N = 0;
            int total_rank = 0;

            List<Integer> ranklist = new ArrayList<>(rankCountMap.values());
            Collections.sort(ranklist, Collections.reverseOrder());


            for (Integer rank : ranklist) {
                int times = rankCountMap.get(rank);
                N += times;
                total_rank += rank*times;
                String perc = String.format(" (%.2f%%)",100.0*times/rankCountMap.size());
                writer.write("rank " + rank + ": " + times + " cases " + perc + "\n");
            }
            double avg_rank = (double)total_rank/N;
            writer.write("average rank " + avg_rank +  ".\n" );
            for (String line : resultlist) {
                writer.write(line + "\n");
            }
            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
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
