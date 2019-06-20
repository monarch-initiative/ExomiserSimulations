#!/usr/bin/env bash

########################    RUN THREES EVALUATOR    #####################
#
# `EselatorSimulations --threes-evaluator` command accepts path to a single phenopacket in the command line (`--pp` option).
# This script runs the

JAR_PATH="simulations-cli/target/simulations-cli-0.2.0-SNAPSHOT.jar"

### --------------------------- VARIABLES --------------------------- ###

# phenopacket directory
PP_DIR=""

# path to spring config location (application.properties)
CONFIG_LOC=""

# path to directory where results will be stored
OUT_DIR=""

# path to VCF file that will be used for simulation of an exome
TEMPLATE_VCF=""

### ----------------------------------------------------------------- ###

check () {
if [[ "${PP_DIR}" = "" ]]; then
  echo "Missing '--pp-dir' argument!"
  usage;
  exit 1
fi

if [[ "${CONFIG_LOC}" = "" ]]; then
  echo "Missing '--config-loc' argument!"
  usage;
  exit 1
fi

if [[ "${OUT_DIR}" = "" ]]; then
  echo "Missing '--out-dir' argument!"
  usage;
  exit 1
fi

if [[ "${TEMPLATE_VCF}" = "" ]]; then
  echo "Missing '--template-vcf' argument!"
  usage;
  exit 1
  fi
}

run () {
# all the phenopacket JSON files
PP_FILES=$(ls $PP_DIR/*.json)

# create command line arguments - `--pp=/path/to/first_phenopacket.json --pp=/path/to/second_phenopacket.json ...`
PP_CLI=""
for F in ${PP_FILES}; do
  PP_CLI="${PP_CLI} --pp=${F}"
done

# now we have everything to run the app
COMMAND="java -jar ${JAR_PATH} --threes-evaluator --spring.config.location=${CONFIG_LOC} --output=${OUT_DIR} --vcf=${TEMPLATE_VCF} ${PP_CLI}"

echo "Running '${COMMAND}'"
$COMMAND
}


usage () {
printf "$(basename ${0%.sh})

USAGE:
    --config-loc       path to directory with Exomiser resources
    --pp-dir           path to directory with phenopackets in JSON format
    --template-vcf     path to VCF file used as a template
    --out-dir          where to write the results

    -h | --help        print this message\n\n"
}


######################## MAIN ####################### MAIN ####################
if [[ "$1" = "" ]]; then
  printf "No parameter inserted!\n"; usage; exit 1
else
  while [[ "$1" != "" ]]; do
    case $1 in
      --config-loc ) shift; CONFIG_LOC=$1;;
      --pp-dir ) shift; PP_DIR=$1;;
      --out-dir ) shift; OUT_DIR=$1;;
      --template-vcf ) shift; TEMPLATE_VCF=$1;;
      -h | --help ) usage; exit 0;;
      *) printf "Unknown parameter: \"$1\"\n"; usage; exit 1;;
    esac
    shift
  done
fi

check;
run;