# Sentieon DNAseq and DNAscope Benchmarking on AWS

This repository contains scripts used for benchmarking Sentieon DNAseq and DNAscope pipelines on AWS EC2.

### Instance Setup

The script at `misc/instance_setup.sh` performs initial setup of the instance and download/installation of software packages used in the benchmark.

### Input datasets

Input datasets for the benchmark are recorded in the `config/config.yaml` (and `config/config_arm.yaml`) files. With the exception of the Element dataset, the pipeline will automatic download input files. Visiting the Element web page recorded in the config file using an internet browser is required for downloading the Element dataset.

### Running benchmarks

The script at `misc/run_benchmarks.sh` was used to run the benchmarks. Before running the scripts, the `SENTIEON_LICENSE` variable inside the script should be properly set.
