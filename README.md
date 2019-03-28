# Eselator simulations

This repo contains code for running simulations using Exomiser version enriched with code for splicing analysis.

## Set up
You need to follow a few steps before you can run the app:
- get and unzip the resource archive - ~35GB ZIP archive containing genome, phenotype, and spliding databases, fasta files, ... 
- update `exomiser.data-directory` inside the *splicing-exomiser-application.properties* file to point to the resource directory

## Run the app

### IntelliJ
The easiest way to run the app is to run it within IntelliJ.
- find class `org.monarchinitiative.eselator.simulations.Play` and add it into the run configurations by clicking on the green triangle next to the class name declaration
- open the newly created configuration and add following into the *Program arguments* field: `--spring.config.location=/home/ielis/soft/exomiser-data/splicing-exomiser-application.properties`
  > the path should point to location of the *splicing-exomiser-application.properties* in the resource directory

After these steps you should be all set.

### Command line
Provide path to phenopacket, VCF file with variants, and file prefix where the results will be written:

```bash
java -jar simulations-cli-0.1.0.jar --pp=example.phenopacket --vcf=example.vcf --output=output/file/prefix
```
