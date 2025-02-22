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

### Simulation mode

In the simulation mode, segmeter generates of intervals (reference) and their corresponding basic and complex queries. This can be used as follows:

```
segmeter sim -o DATADIR [-h] [-n INVLNUMS] [-m MAX_CHROMLEN] [-c SIMNAME] [-g GAPSIZE] [-i INTVLSIZE]

```

| Argument | Description |
| -------- | ----------- |
| -o, --datadir | output folder for the benchmark/simulation results. Note this also serves as input folder for the benchmarking |
| -n, --intvlnums | Number of intervals to simulate (should be divisible by 10). Can be a comma separated list of intervals (for different datasets). Can be abbreviated for thousands, millions (e.g., 10K, 1M). Default is 10.|
| -m, --max_chromlen | maximum length (in base pairs) of the simulated chromosomes. The default maximum length is set to one billion (e.g., 1000000000). In this is exceeded in the simulation, segmeter creates new scaffolds. |
| -c, --simname | name of the simulation, used for the output folder |
| -g, --gapsize | random size of the gaps (min and max) between the intervals. Default is 100-5000 |
| -i, --intvlsize | random size (min and max) of the intervals. Default is 100-10000 |

This will generates output files in the 'DATADIR/simname/BED' folder. The files are in BED format (currently the only supported format) and can be used for benchmarking. In particular, the filers are located in the following folders:

```
DATADIR/simname/BED/ref/ # reference intervals
DATADIR/simname/BED/basic/ # basic queries
DATADIR/simname/BED/complex/ # complex queries
```

#### Reference

`DATADIR/simname/BED/ref` contains the reference intervals. This is a BED4 file with the interval ID in the fourth column.
For each value in `INTVLNUMS`, there is a corresponding file with the intervals.

#### Basic queries

`DATADIR/simname/BED/basic` contains the basic queries. In total, segmeter generates ten basic queries for each interval in the reference.
For each value in `INTVLNUMS`, there is a corresponding folder with the basic queries. In each folder, there is a  subfolders for the different
types of basic queries:
```
DATADIR/simname/BED/basic/perfect/ # perfect overlaps (boundaries are identical with reference)
DATADIR/simname/BED/basic/5p-partial/ # partial overlaps on 5' end
DATADIR/simname/BED/basic/3p-partial/ # partial overlaps on 3' end
DATADIR/simname/BED/basic/contained/ # overlap is contained within the reference
DATADIR/simname/BED/basic/enclosed/ # overlap encloses the reference
DATADIR/simname/BED/basic/perfect-gap/ # perfect overlap with a gap between reference intervals
DATADIR/simname/BED/basic/left-adjacent-gap/ # adjacent to the 5'-end of the interval (no overlap)
DATADIR/simname/BED/basic/right-adjacent-gap/ # adjacent to the 3'-end of the interval (no overlap)
DATADIR/simname/BED/basic/mid-gap1/ # random overlap with a gap
DATADIR/simname/BED/basic/mid-gap2/ # random overlap with a gap
```

In each of the subfolders, there is a BED4 file for each of the specified INTVLNUMS (e.g., `DATADIR/simname/BED/basic/query/<query_type>/<INTVLNUM>.bed`).
In addition, the queries are subsample to 10-100% of the queries and stored in corresponding files (e.g., `DATADIR/simname/BED/basic/query/<query_type>/INTVLNUM_<PERCENT>p.bed`).

##### Truth

In addition, the truth files are stored in `DATADIR/simname/BED/basic/truth/<INTVLNUM>.bed`. These files contain the queries and their corresponding reference intervals. In each line
the query interval is followed by the reference interval. The reference interval is the interval that the query should overlap with. In the fourth column, the combined ID of the query and
reference interval is stored.

```
chr13	309	1623	chr13	584	4573	intvl_1_5p:intvl_1
chr13	3221	4662	chr13	584	4573	intvl_1_3p:intvl_1
chr13	1823	2593	chr13	584	4573	intvl_1_contained:intvl_1
chr13	519	4649	chr13	584	4573	intvl_1_enclosed:intvl_1
chr13	584	4573	chr13	584	4573	intvl_1_perfect:intvl_1
```

#### Complex queries

`DATADIR/simname/BED/complex` contains the complex queries.









### Benchmark mode

In the benchmark mode, segmeter reads a dataset of intervals from a file and evaluates the performance of a given interval retrieval algorithm. This can be used as follows:






```
usage: main.py [-h] -o DATADIR [-f {BED}] [-n INTVLNUMS] [-s SUBSET]
               [-m MAX_CHROMLEN] [-b BENCHNAME] [-c SIMNAME]
               [-t {tabix,bedtools,bedtools_sorted,bedtools_tabix,bedops,bedmaps,giggle,granges,gia,gia_sorted,bedtk,bedtk_sorted,igd,ailist,ucsc}]
               [-g GAPSIZE] [-i INTVLSIZE]
               {sim,bench}
```


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

In additon, we provide a ready-to-use Docker container that has segmeter preconfigured. It can be found at [dockerhub](https://hub.docker.com/r/yanglabinfo/segmeter). We provide three different containers that can be used for the different tools.

| Container      | Tools      | Container tag |
| ------------- | ------------- | ------------- |
| giggle | giggle | segmeter:giggle-latest |
| others | AIList, BEDops, bedtools, bedtk, IGD, tabix, UCSC utils | segmeter:others-latest |
| rust-tools | gia, granges | segmeter:rust-tools-latest |

This can used with the following commands:
```
docker run -it -d -v /folder/on/host/:/folder/in/container/ yanglabinfo/segmeter:<container_tag> /bin/bash
docker exec <container_id> segmeter <args>
```


## Singularity

Segmeter is provided as docker container that can be used using `docker pull yanglabinfo/segmeter`. Consequently, this can be also used
