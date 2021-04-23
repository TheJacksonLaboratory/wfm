#!/bin/bash

conda env create -f environment.yml
wget -P cromwell/ https://github.com/broadinstitute/cromwell/releases/download/54/cromwell-54.jar
cd nextflow/
wget -qO- https://get.nextflow.io | bash
cd ../
sed -i "s,CURRENT,$(pwd),g" cromwell/inputs.json
sed -i "s,CURRENT,$(pwd),g" nextflow/nextflow.config
sed -i "s,CURRENT,$(pwd),g" snakemake/config.yaml
