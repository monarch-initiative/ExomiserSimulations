# Eselator simulations

This repo contains code for running simulations using Exomiser version enriched with code for splicing analysis. The code will receive a `Phenopacket`, simulate an exome VCF and run the Exomiser analysis.

Simulated exome VCF file is being created in a naive way at the moment. A single VCF file with variants is required and the variants from `Phenopacket` are spiked in between the variants.

## Set up
You need to follow a few steps before you can run the app:
- get and unzip the resource archive - ~35GB ZIP archive containing genome, phenotype, and spliding databases, fasta files, ... 
- update `exomiser.data-directory` inside the *splicing-exomiser-application.properties* file to point to the resource directory

## Inject phenopacket variants & phenotype into VCF file and run Exomiser analysis 

### Command line
Provide path to phenopacket, VCF file with variants, and file prefix where the results will be written:

```bash
java -jar simulations-cli-0.1.0.jar --exomiser.data-directory=/path/to/exomiser-data/directory --pp=/path/to/phenopacket.json --vcf=/path/to/vcf --output=/path/to/output
```
> **Note:** 
  - this repo uses Exomiser with splicing code at the moment, so `--exomiser.data-directory` must point to exomiser data bundle where the splicing files are present as well.
  - `--pp` must point to a `Phenopacket` containing a single individual. Familial analysis is not supported at the moment. Only **present** (not negated) `Phenotype`s are injected into Exomiser analysis.
  - `--output` is supposed to be a file prefix, result files will be created for all supported inheritance modes and output types (TSV, VCF, HTML, ...)

### IntelliJ
The above analysis can also be run within IntelliJ. You just have to create a run configuration:
- find class `org.monarchinitiative.eselator.simulations.Play` and add it into the run configurations by clicking on the green triangle next to the class name
- open the newly created configuration and add the following into *Program arguments* field: `--exomiser.data-directory=/path/to/exomiser-data/directory --pp=/path/to/phenopacket.json --vcf=/path/to/vcf --output=/path/to/output`
- hit the *Run* button and hope for the best :)

## Limitations
There are some hardcoded constraints present at the moment:

- only runs with *RefSeq* transcripts and *hg19* genome assembly
- only single sample analysis is supported
- only one *naive* exome simulation method is available at the moment
- exomiser version with splicing code (`11.0.0-SP-1`) must be installed in your local Maven repo in order to compile & run the app
