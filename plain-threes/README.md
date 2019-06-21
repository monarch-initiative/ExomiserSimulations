# Plain threes

This module performs simulations/evaluation of 3S code.

**Available commands:**
- [Score phenopackets](#Score-phenopackets) - apply all 3S scoring strategies to score variants in given phenopackets and write results into a TSV file
- [Simulate case and run Exomiser](#Simulate-case-and-run-Exomiser) - take a directory of Phenopacket and simulate exome VCF for each one. Then run Exomiser either with or without SPLICING score. Store ranks of causal genes in TSV file and save Exomiser results (HTML, TSV, etc..)
- [Move phenopackets without phenotype](#Move-phenopackets-without-phenotype) - some phenopackets contain 0 HPO terms which will crash Exomiser analysis where we use HiPhive prioritiser. This command will move such Phenopackets into separate directory.


## Run all
This command runs `--simulate-case-and-run-exomiser` and `--score-phenopackets`.


```bash
java -jar plain-threes-0.2.0-SNAPSHOT.jar
--simulate-case-and-run-exomiser
--score-phenopackets
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets/dir
--vcf=/path/to/template.vcf
--output-exomiser=/path/to/output-dir
--output-scores=/path/to/results.tsv
```

## Score phenopackets

```bash
java -jar plain-threes-0.2.0-SNAPSHOT.jar
--score-phenopackets
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets/dir
--output-scores=/path/to/results.tsv
``` 
> Note: You can also specify path to individual phenopackets using `--pp` option.
 

## Simulate case and run Exomiser

```bash
java -jar plain-threes-0.2.0-SNAPSHOT.jar
--simulate-case-and-run-exomiser
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets/dir
--vcf=/path/to/template.vcf
--output-exomiser=/path/to/output-dir
```
> Note: You can also specify path to individual phenopackets using `--pp` option.


## Move phenopackets without phenotype

```bash
java -jar plain-threes-0.2.0-SNAPSHOT.jar
--move-phenopackets-without-phenotype
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets-dir
```
