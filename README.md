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

In addition, segmeter generates for each specified `INTVLNUM`, a file with the length of each simulated chromosome (`DATADIR/simname/BED/<INTVLNUM>_chrlens.txt`),
and the number of intervals per chromosome (`DATADIR/simname/BED/<INTVLNUM>_chrnums.txt`).

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

`DATADIR/simname/BED/complex` contains the complex queries. Currently, this only includes `mult` queries which basically cover multiple reference intervals. According to the number of intervals that
are covered in a complex query, the queries are stored in deciles (e.g., `DATADIR/simname/BED/complex/query/mult/<INTVLNUM>_<DECILE>bin.bed`). Again this contains the queries in BED4 format with and
identifier in the fourth column. Note that `mult_13` indicates that this query covers 13 reference intervals:
```
chr12	45837	160962	mult_13
chr12	146101	264093	mult_14
chr12	87719	214921	mult_15
chr12	122976	259591	mult_16
chr14	20381	105691	mult_13
```

##### Truth

The truth files for the complex queries are stored in `DATADIR/simname/BED/complex/truth/<INTVLNUM>.bed`. This consists of the query interval and the corresponding number of intervals that are covered by the query.
```
chr11	43495	63230	mult_2	2
chr11	106913	119696	mult_3	3
chr11	90920	113986	mult_4	4
chr11	136183	172964	mult_5	5
chr11	236310	283001	mult_6	6
```

### Benchmark mode

In the benchmark mode, segmeter reads a dataset of intervals from a file and evaluates the performance of a given interval retrieval algorithm. This can be used as follows:
```
segmeter bench -o DATADIR [-h] [-n INTVLNUMS] [-s SUBSET] [-b BENCHNAME] [-c SIMNAME] [-t TOOL]
```

| Argument | Description |
| -------- | ----------- |
| -o, --datadir | input/output folder for the simulation results. Note that this folder must contains a subfolder `sim` that contains the simulated interval data |
| -n, --intvlnums | Number of intervals to benchmark. When multiple datasets are benchmark, this should be a comma separated list (same in in simulation). Note that this should have been simulated before. |
| -s, --subset | subset (in percentage) of the intervals to use for benchmarking. Format should be either XX-YY or XX,YY-ZZ. If this is left empty, all subsets/deciles are used |
| -b, --benchname | name of the benchmark, used for the output folder. This allows to perform multiple benchmarks |
| -c, --simname | name of the simulation data that is being used. Note that this should be the same as the name of the simulation data that was used for the simulation |
| -t, --tool | tool to benchmark. Currently, the following tools are supported: `tabix`, `bedtools`, `bedtools_sorted`, `bedtools_tabix`, `bedops`, `bedmaps`, `giggle`, `granges`, `gia`, `bedtk`, `bedtk_sorted`, `igd`, `ailist`, `ucsc` |

Note that `bedtools_sorted` and `bedtk_sorted` are the same as `bedtools` and `bedtk`, respectively, but the input files are sorted before the benchmarking.
In the case of `bedtools_tabix`, the input files are sorted and indexed using `tabix` for random access and the queried using `bedtools`.

This generates a separate output folder for each benchmark tool in the folder `DATADIR/bench/benchname/` with a subfolder for each INTVLNUM.
In additional subfolders (`precision` and `stats`), the precision and statistics are stored.
The precision is stored in a file `DATADIR/bench/benchname/precision/<INTVLNUM>_<PERCENT>.txt` and the
statistics in `DATADIR/bench/benchname/stats/<INTVLNUM>_<PERCENT>.txt`.

#### Precision

In the precision files, the precision of the tool on basic and complex queries is stored separately:
```
intvlnum	subset	TP	FP	TN	FN	Precision	Recall	F1
1000	10%	500	0	500	0	1.0	1.0	1.0

intvlnum	bin	distance
1000	10bin	0
```

The upper part of the file contains the precision, recall, and F1 score for the basic queries and subset (e.g., 10% of the queries).
The lower part contains the distance which is the absolute difference between expected and observed number of intervals covered by the complex query.
Note this only represents a decile (e.g., 10bin), in other words, the queries that cover 10% of the reference intervals per chromosome.

#### Statistics

In the statistics files, the statistics of the tool on basic and complex queries is stored separately:
```
intvlnum	data_type	query_type	time	max_RSS(MB)
1000	basic	perfect_100%	0.00279	1.2265625
1000	basic	5p-partial_100%	0.00261	1.22265625
1000	basic	3p-partial_100%	0.00282	1.1640625
1000	basic	enclosed_100%	0.00255	1.2265625
1000	basic	contained_100%	0.00294	1.22265625
1000	basic	perfect-gap_100%	0.00329	1.2265625
1000	basic	left-adjacent-gap_100%	0.00244	1.1640625
1000	basic	right-adjacent-gap_100%	0.00211	1.1640625
1000	basic	mid-gap1_100%	0.00215	1.2265625
1000	basic	mid-gap2_100%	0.0022	1.2265625
1000	complex	mult_100bin	0.00193	1.2265625
```

The file contains the time and memory usage of the tool for each query type and subset. The time is in seconds and the memory usage in MB.

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
