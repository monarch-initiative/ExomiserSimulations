# Plain threes

This module performs simulations/evaluation of 3S code.

**Available commands:**
- [Score phenopackets](#Score-phenopackets) - apply all 3S scoring strategies to score variants in given phenopackets and write results into a TSV file
- [Simulate case and run Exomiser](#Simulate-case-and-run-Exomiser) - take a directory of Phenopacket and simulate exome VCF for each one. Then run Exomiser either with or without SPLICING score. Store ranks of causal genes in TSV file and save Exomiser results (HTML, TSV, etc..)
- [Clinvar scorer](#Clinvar-scorer) - select variants with benign or likely benign clinical significance (see `--strict` flag) and score variants using all splicing strategies. Write the results into a TSV file
- [Move phenopackets without phenotype](#Move-phenopackets-without-phenotype) - some phenopackets contain 0 HPO terms which will crash Exomiser analysis where we use HiPhive prioritiser. This command will move such Phenopackets into separate directory

## Run all
This command runs `--simulate-case-and-run-exomiser` and `--score-phenopackets`.

**The point** of running both command is that you can evaluate *n* selected cases by Exomiser and also see how the individual splicing scorers did perform for each splicing variant.

This way Exomiser result files will be created for each phenopacket in the `--output-exomiser` directory. Scores produced by each scorer will be written into TSV file at `--output-scores` path.

```bash
java -jar plain-threes-0.2.1.jar
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
java -jar plain-threes-0.2.1.jar
--score-phenopackets
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets/dir
--output-scores=/path/to/results.tsv
```
> Note: You can also specify path to individual phenopackets using `--pp` option.


## Simulate case and run Exomiser

```bash
java -jar plain-threes-0.2.1.jar
--simulate-case-and-run-exomiser
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets/dir
--vcf=/path/to/template.vcf
--output-exomiser=/path/to/output-dir
```
> Note: You can also specify path to individual phenopackets using `--pp` option.

## Clinvar scorer

```bash
java -jar plain-threes-0.2.1.jar 
--clinvar-scorer
--spring.config.location=/path/to/application.properties
--clinvar-vcf=/path/to/clinvar.vcf.gz
--output-clinvar=/path/to/output_file.tsv
--strict # if you want to only process the benign variants, not likely benign
```
> Note: I did not test this command but is should run
**Limitations:**

- Does not work with variants from other chromosomes than \[1..22,X,Y\] (e.g. `MT`) at the moment

## Move phenopackets without phenotype

```bash
java -jar plain-threes-0.2.1.jar
--move-phenopackets-without-phenotype
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets-dir
```
