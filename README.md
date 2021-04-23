# Introduction to Scientific Workflow Managers

This repository provides a starting point to using three popular scientific workflow managers - Nextflow, Cromwell, and Snakemake. This material is used as instructional material for the Research IT workshop. 

## Getting Started

In order to use these workflow managers, first download Miniconda3. This can be done on the command line (assuming you are using our Linux cluster) by doing the following. First, navigate to a directory in which you have at least a few GB of available space. Then, use the following commands.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Once Miniconda is installed, you can then download the required software. Clone the repository, and then run `install.sh` to install the pre-made conda environment.

```
git clone https://github.com/TheJacksonLaboratory/wfm.git
cd wfm
bash install.sh
```

This will create a conda environment named `wfm`, as well as download the Nextflow and Cromwell execution engines into their designated folders. The last step will be 

