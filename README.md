# Exomiser simulations

This repo contains code for running simulations using Exomiser version enriched with code for splicing analysis. The code will receive a `Phenopacket`, simulate an exome VCF and run the Exomiser analysis.

Simulated exome VCF file is being created in a naive way at the moment. A single VCF file with variants is required and the variants from `Phenopacket` are spiked in between the variants.

## Inject phenopacket variants & phenotype into VCF file and run Exomiser analysis 

**Command line arguments:**
- `--single-vcf-simulation` - flag required to run the command
- `--exomiser.data-directory` - path to Exomiser data bundle
- `--pp` - path to Phenopacket in JSON format
- `--vcf` - path to VCF file with variants, where variants from Phenopacket will be spiked in
- `--output` - prefix of the Exomiser result files

Example run:
```bash
java -jar simulations-cli-0.1.0.jar --single-vcf-simulation --exomiser.data-directory=/path/to/exomiser-data/directory --pp=/path/to/phenopacket.json --vcf=/path/to/vcf --output=/path/to/output
```

**Limitations:**
There are some hardcoded constraints present at the moment:

- only runs with *RefSeq* transcripts and *hg19* genome assembly
- only single sample analysis is supported
- only one *naive* exome simulation method is available at the moment
- Exomiser version with splicing code (`11.0.0-SP-1`) must be installed in your local Maven repo in order to compile & run the app
- no information regarding expected mode of inheritance is extracted from the `Phenopacket` 

## Run LIRICAL

**Command line arguments:**
- `--lirical` - flag required to run the command
- `--exomiser.data-directory` - path to Exomiser data bundle
- `--pp` - path to directory with Phenopackets in JSON format. Note: phenopackets **only** must be present in the directory
- `--vcf` - path to VCF file with variants, where variants from Phenopacket will be spiked in
