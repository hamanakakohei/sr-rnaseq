#!/bin/bash
set -euo pipefail

WDL=workflows/03_rnaseq_pipeline_fastq.edit.wdl
OPTIONS=cromwell-options.json

for INPUT in inputs/03/*.json; do
    SAMPLE=$(basename $INPUT .json)
    echo "Running $SAMPLE ..."
    java -jar cromwell.jar run $WDL -i $INPUT -o $OPTIONS
done


