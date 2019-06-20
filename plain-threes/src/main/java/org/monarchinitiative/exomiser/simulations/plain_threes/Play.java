package org.monarchinitiative.exomiser.simulations.plain_threes;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

/**
 *
 */
@SpringBootApplication
public class Play {


    public static void main(String[] args) {
        SpringApplication.run(Play.class, args);
    }

/*
    @Bean
    public SplicingEvaluator splicingEvaluator(SplicingTranscriptLocator splicingTranscriptLocator, ScorerFactory scorerFactory) {
        return new SimpleSplicingEvaluator(splicingTranscriptLocator, scorerFactory);
    }

    @Bean
    public SplicingTranscriptLocator splicingTranscriptLocator(SplicingParameters splicingParameters, GenomeCoordinatesFlipper genomeCoordinatesFlipper) {
        return new NaiveSplicingTranscriptLocator(splicingParameters, genomeCoordinatesFlipper);
    }

    @Bean
    public SplicingParameters splicingParameters(SplicingInformationContentAnnotator splicingInformationContentAnnotator) {
        return splicingInformationContentAnnotator.getSplicingParameters();
    }

    @Bean
    public GenomeCoordinatesFlipper genomeCoordinatesFlipper(Map<String, Integer> contigLengthMap) {
        return new GenomeCoordinatesFlipper(contigLengthMap);
    }

*/


}
