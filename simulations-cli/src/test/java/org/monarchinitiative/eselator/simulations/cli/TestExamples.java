package org.monarchinitiative.eselator.simulations.cli;

import com.google.protobuf.Timestamp;
import org.phenopackets.schema.v1.core.MetaData;
import org.phenopackets.schema.v1.core.OntologyClass;
import org.phenopackets.schema.v1.core.Resource;

import java.time.Instant;
import java.util.Arrays;
import java.util.stream.Collectors;

public class TestExamples {

    public static final OntologyClass HOMO_SAPIENS = ontologyClass("NCBITaxon:9606", "Homo sapiens");

    public static final OntologyClass TAS = ontologyClass("ECO:0000033", "author statement supported by traceable reference");

    /**
     * Source https://bioportal.bioontology.org/ontologies/GENO/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGENO_0000135
     */
    public static final OntologyClass HET = ontologyClass("GENO:0000135", "heterozygous");

    /**
     * Source https://bioportal.bioontology.org/ontologies/GENO/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGENO_0000136
     */
    public static final OntologyClass HOM_ALT = ontologyClass("GENO:0000136", "homozygous");

    public static final OntologyClass HEMIZYGOUS = ontologyClass("GENO:0000134", "hemizygous");


    public static OntologyClass ontologyClass(String id, String label) {
        return OntologyClass.newBuilder()
                .setId(id)
                .setLabel(label)
                .build();
    }


    public static MetaData getMetaData(Resource... resources) {
        return MetaData.newBuilder()
                .setCreated(Timestamp.newBuilder()
                        .setSeconds(Instant.parse("2019-03-20T18:14:54Z").getEpochSecond())
                        .build())
                .setCreatedBy("Mr Fantastic")
                .setSubmittedBy("Test corp.")
                .addAllResources(Arrays.stream(resources).collect(Collectors.toList()))
                .build();
    }

    public static Resource getNcbiOrganismalClassificationResource() {
        return Resource.newBuilder()
                .setId("ncbitaxon")
                .setName("NCBI organismal classification")
                .setNamespacePrefix("NCBITaxon")
                .setIriPrefix("http://purl.obolibrary.org/obo/NCBITaxon_")
                .setUrl("http://purl.obolibrary.org/obo/ncbitaxon.owl")
                .setVersion("19-03-2019")
                .build();
    }

    public static Resource getEvidenceOntologyResource() {
        return Resource.newBuilder()
                .setId("eco")
                .setName("Evidence ontology")
                .setNamespacePrefix("ECO")
                .setIriPrefix("http://purl.obolibrary.org/obo/ECO_")
                .setUrl("http://purl.obolibrary.org/obo/eco.owl")
                .setVersion("19-03-2019")
                .build();
    }

    public static Resource getGenotypeOntologyResource() {
        return Resource.newBuilder()
                .setId("geno")
                .setName("Genotype Ontology")
                .setNamespacePrefix("GENO")
                .setIriPrefix("http://purl.obolibrary.org/obo/GENO_")
                .setUrl("http://purl.obolibrary.org/obo/geno.owl")
                .setVersion("19-03-2018")
                .build();
    }

    public static Resource getHPOResource() {
        return Resource.newBuilder()
                .setId("hp")
                .setName("human phenotype ontology")
                .setNamespacePrefix("HP")
                .setIriPrefix("http://purl.obolibrary.org/obo/HP_")
                .setUrl("http://purl.obolibrary.org/obo/hp.owl")
                .setVersion("2018-03-08")
                .build();
    }
}