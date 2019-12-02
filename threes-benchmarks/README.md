# 3S benchmarks

This module performs simulations/evaluation of 3S code alone.

**Available commands:**
- [Score phenopackets](#Score-phenopackets) - apply all 3S scoring strategies to score variants in given phenopackets and write results into a TSV file
- [Clinvar scorer](#Clinvar-scorer) - select variants with benign or likely benign clinical significance (see `--strict` flag) and score variants using all splicing strategies. Write the results into a TSV file


## Score phenopackets

```bash
java -jar plain-threes-0.2.1.jar
--score-phenopackets
--spring.config.location=/path/to/application.properties
--pp-dir=/path/to/phenopackets/dir
--output-scores=/path/to/results.tsv
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

