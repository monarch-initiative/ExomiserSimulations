# Plain threes

This module performs simulations/evaluation of 3S code when used as plugin for Exomiser.

**Available commands:**
- [Simulate case and run Exomiser](#Simulate-case-and-run-Exomiser) - take a directory of Phenopacket and simulate exome VCF for each one. Then run Exomiser either with or without SPLICING score. Store ranks of causal genes in TSV file and save Exomiser results (HTML, TSV, etc..)
- [Move phenopackets without phenotype](#Move-phenopackets-without-phenotype) - some phenopackets contain 0 HPO terms which will crash Exomiser analysis where we use HiPhive prioritiser. This command will move such Phenopackets into separate directory

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

## Move phenopackets without phenotype

```bash
java -jar plain-threes-0.2.1.jar
--move-phenopackets-without-phenotype
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets-dir
```


---------

## Run all (not working anymore)
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
