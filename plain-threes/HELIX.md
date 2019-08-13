# Helix

This is how TORQUE batch script looks like:

```bash
#!/bin/bash

#PBS -N 3S_simulations
#PBS -l nodes=1,mem=24gb,walltime=16:30:00,file=20gb
#PBS -q batch
#PBS -V
#PBS -M daniel.danis@jax.org
#PBS -m abe

module load java/1.8.0_73


RUN_ID="20190812"
ES_JAR_PATH="/home/danisd/bin/jars/plain-threes-0.2.1.jar"
ES_APP_PROP="/projects/robinson-lab/3S/exomiser-data/helix-splicing-exomiser-application.properties"
PP_DIR_PATH="/projects/robinson-lab/3S/${RUN_ID}-phenopackets-checked"
GIAB_VCF="/projects/robinson-lab/3S/GIAB/project.NIST.hc.snps.indels.NIST7035.vcf"
OUTPUT_EXOMISER="/projects/robinson-lab/3S/${RUN_ID}-output-checked"
OUTPUT_PP="${OUTPUT_EXOMISER}/phenopacket_scores.tsv"


cd ${PBS_O_WORKDIR}

java -jar ${ES_JAR_PATH} --simulate-case-and-run-exomiser --score-phenopackets --spring.config.location=${ES_APP_PROP} --pp-dir=${PP_DIR_PATH} --vcf=${GIAB_VCF} --output-exomiser=${OUTPUT_EXOMISER} --output-scores=${OUTPUT_PP}
```

