<div align="left">
    <h1>segmeter</h1>
    <img src="https://img.shields.io/github/v/release/ylab-hi/segmeter">
    <img src="https://github.com/ylab-hi/ScanNeo2/actions/workflows/linting.yml/badge.svg" alt="Workflow status badge">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg">
    <img src="https://img.shields.io/github/downloads/ylab-hi/segmeter/total.svg">
    <img src="https://img.shields.io/github/contributors/ylab-hi/segmeter">
    <img src="https://img.shields.io/github/last-commit/ylab-hi/segmeter">
    <img src="https://img.shields.io/github/commits-since/ylab-hi/segmeter/latest">
    <img src="https://img.shields.io/github/stars/ylab-hi/segmeter?style=social">
    <img src="https://img.shields.io/github/forks/ylab-hi/segmeter?style=social">
</div>

## What is segmeter

This is a tool for simulating interval data and benchmarking tool for interval retrieval.

## Usage

segmeter currently supports two modes of operation: `sim` (e.g., simulate) and `bench` (e.g., benchmark).
In the `sim` mode, segmeter generates a synthetic dataset of intervals and writes it to a file. In the `bench` mode,
segmeter reads a dataset of intervals from a file and evaluates the performance of a given interval retrieval algorithm.

usage: main.py [-h] -o DATADIR [-f {BED}] [-n INTVLNUMS] [-s SUBSET]
               [-m MAX_CHROMLEN] [-b BENCHNAME] [-c SIMNAME]
               [-t {tabix,bedtools,bedtools_sorted,bedtools_tabix,bedops,bedmaps,giggle,granges,gia,gia_sorted,bedtk,bedtk_sorted,igd,ailist,ucsc}]
               [-g GAPSIZE] [-i INTVLSIZE]
               {sim,bench}

Benchmarking tool for interval files

positional arguments:
  {sim,bench}           modus in benchmarking

options:
  -h, --help            show this help message and exit
  -o DATADIR, --datadir DATADIR
                        output folder for the benchmark/simulation results.
                        Note this also serves as input folder for the
                        benchmarking
  -f {BED}, --format {BED}
                        format of the files to benchmark
  -n INTVLNUMS, --intvlnums INTVLNUMS
                        Number of intervals to simulate (should be divisible
                        by 10).
  -s SUBSET, --subset SUBSET
                        subset of the intervals to use for benchmarking.
                        Format should be either XX-YY or XX,YY-ZZ
  -m MAX_CHROMLEN, --max_chromlen MAX_CHROMLEN
                        maximum length of the simulated chromosomes
  -b BENCHNAME, --benchname BENCHNAME
                        name of the benchmark, used for the output folder
  -c SIMNAME, --simname SIMNAME
                        name of the simulation, used for the output folder.
                        Only used in simulation mode
  -t {tabix,bedtools,bedtools_sorted,bedtools_tabix,bedops,bedmaps,giggle,granges,gia,gia_sorted,bedtk,bedtk_sorted,igd,ailist,ucsc}, --tool {tabix,bedtools,bedtools_sorted,bedtools_tabix,bedops,bedmaps,giggle,granges,gia,gia_sorted,bedtk,bedtk_sorted,igd,ailist,ucsc}
                        tool to benchmark
  -g GAPSIZE, --gapsize GAPSIZE
                        random size of the gaps (min and max) between the
                        intervals
  -i INTVLSIZE, --intvlsize INTVLSIZE
                        random size (min and max) of the intervals



## Docker

In additon, we provide a ready-to-use Docker container that has segmeter preconfigured. It can be found at [dockerhub](https://hub.docker.com/r/yanglabinfo/segmeter)



Due to compability issues, the tools are available in different containers



## Singularity

Segmeter is provided as docker container that can be used using `docker pull yanglabinfo/segmeter`.
