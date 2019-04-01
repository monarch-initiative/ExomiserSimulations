#!/usr/bin/env bash

########################    RUN SPLICE SCORER    ########################
#
# `EselatorSimulations --splice-scorer` command accepts path to a single phenopacket in the command line (`--pp` option).
# This script runs the

JAR_PATH="simulations-cli/target/simulations-cli-0.2.0-SNAPSHOT.jar"

### --------------------------- VARIABLES --------------------------- ###

# phenopacket directory
PP_DIR=""

# exomiser data directory
EXOMISER_DATADIR=""

# output file path
OUT_FILE=""

# strategy - advanced by default
STRATEGY="advanced"

### ----------------------------------------------------------------- ###

check () {
if [[ "${PP_DIR}" = "" ]]; then
  echo "Missing '--pp-dir' argument!"
  usage;
  exit 1
fi

if [[ "${EXOMISER_DATADIR}" = "" ]]; then
  echo "Missing '--exomiser-dir' argument!"
  usage;
  exit 1
fi

if [[ "${OUT_FILE}" = "" ]]; then
  echo "Missing '--out-file' argument!"
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
COMMAND="java -jar ${JAR_PATH} --splice-scorer --exomiser.data-directory=${EXOMISER_DATADIR} --strategy=${STRATEGY} --output=${OUT_FILE} ${PP_CLI}"

echo "Running '${COMMAND}'"
$COMMAND
}


usage () {
printf "$(basename ${0%.sh})

USAGE:
    --exomiser-dir     path to directory with Exomiser resources
    --pp-dir           path to directory with phenopackets in JSON format
    --out-file         where to write the results
    --strategy         scoring strategy [advanced] {sigmoid,advanced}
    -h | --help        print this message\n\n"
}


######################## MAIN ####################### MAIN ####################
if [[ "$1" = "" ]]; then
  printf "No parameter inserted!\n"; usage; exit 1
else
  while [[ "$1" != "" ]]; do
    case $1 in
      --exomiser-dir ) shift; EXOMISER_DATADIR=$1;;
      --pp-dir ) shift; PP_DIR=$1;;
      --out-file ) shift; OUT_FILE=$1;;
      --strategy ) shift; STRATEGY=$1;;
      -h | --help ) usage; exit 0;;
      *) printf "Unknown parameter: \"$1\"\n"; usage; exit 1;;
    esac
    shift
  done
fi

check;
run;