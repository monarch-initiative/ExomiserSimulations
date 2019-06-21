# Exomiser simulations

This repo contains code for running simulations using Exomiser. The code usually receives one or more Phenopackets, simulate an exome and run the Exomiser analysis.

Simulated exome VCF file is being created in a naive way at the moment. A single VCF file with variants is required and the variants from `Phenopacket` are spiked in between the variants.

## Run *Lirical* simulations

You have to provide following arguments to run the code:

**Command line argumets:**
```bash
java -jar simulations-cli-0.2.0-SNAPSHOT.jar 
--lirical 
--spring.config.location=/path/to/application.properties 
--pp=/path/to/phenopacket-dir 
--vcf=/path/to/template.vcf
```

- `--spring.config.location` should point to *application.properties* file with paths to Exomiser resources

## Run *3S*-related simulations

See the README file in the `plain-threes` module.

## Run splice scorers on ClinVar VCF

**Command line arguments:**
- `--clinvar-scorer` - flag required to run the command
- `--exomiser.data-directory` - path to Exomiser data bundle - REQUIRED
- `--clinvar-vcf` - path to ClinVar VCF file, the file should contain only variants on chromosomes \[1..22,X,Y\] - REQUIRED
- `--output` - path where result TSV file will be written - REQUIRED

Example run:
```bash
java -jar simulations-cli-0.2.0.jar --clinvar-scorer --exomiser.data-directory=/path/to/exomiser-data/directory --clinvar-vcf=/path/to/clinvar_vcf.gz  --output=/path/to/output_file.tsv
```

**Limitations:**

- Does not work with variants from other chromosomes than \[1..22,X,Y\] (e.g. `MT`) at the moment   
