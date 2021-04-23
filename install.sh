#!/bin/bash

conda env create -f environment.yml
wget -P cromwell/ https://github.com/broadinstitute/cromwell/releases/download/54/cromwell-54.jar
cd nextflow/
wget -qO- https://get.nextflow.io | bash
cd ../
