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
segmeter sim -o DATADIR [-h] [-n INVLNUMS] [-m MAX_CHROMLEN] [-c SIMNAME] [-g GAPSIZE] [-i INTVLSIZE] [--max_span MAX_SPAN]

```

| Argument | Description |
| -------- | ----------- |
| -o, --datadir | output folder for the benchmark/simulation results. Note this also serves as input folder for the benchmarking |
| -n, --intvlnums | Number of intervals to simulate (should be divisible by 10). Can be a comma separated list of intervals (for different datasets). Can be abbreviated for thousands, millions (e.g., 10K, 1M). Default is 10.|
| -m, --max_chromlen | maximum length (in base pairs) of the simulated chromosomes. The default maximum length is set to one billion (e.g., 1000000000). In this is exceeded in the simulation, segmeter creates new scaffolds. |
| -c, --simname | name of the simulation, used for the output folder |
| -g, --gapsize | random size of the gaps (min and max) between the intervals. Default is 100-5000 |
| -i, --intvlsize | random size (min and max) of the intervals. Default is 100-10000 |
| --max_span | maximum number of reference intervals that a complex query covers (a multiple of 10, at least 10, so that every span falls into one of the ten bins). Not limited by default. The output of the complex queries grows quadratically with the number of intervals per chromosome, so a limit (e.g., 1000) is recommended for more than 10K intervals |

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

`DATADIR/simname/BED/complex` contains the complex queries. Currently, this only includes `mult` queries which basically cover multiple reference intervals. For each chromosome, there is one query
for each number of covered intervals, from 2 up to the number of intervals on the chromosome or `--max_span`, whichever is smaller. According to the number of intervals that
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
segmeter bench -o DATADIR -t TOOL [-h] [-r] [-n INTVLNUMS] [-s SUBSET] [-b BENCHNAME] [-c SIMNAME] [--query QUERY] [--target TARGET]
```

| Argument | Description |
| -------- | ----------- |
| -r, --simdata | use simulated data for benchmarking. Note that this will use the simulated data (mode sim) as input for the benchmarking |
| --query | query file used for benchmarking (not used when benchmarking simulated data) |
| --target | target file used for benchmarking (not used when benchmarking simulated data) |
| -o, --datadir | input/output folder. With `-r`, it must contain the subfolder `sim` with the simulated interval data; with `--target`/`--query`, the target is copied to `DATADIR/ref/` |
| -n, --intvlnums | Number of intervals to benchmark. When multiple datasets are benchmark, this should be a comma separated list (same in in simulation). Note that this should have been simulated before. |
| -s, --subset | subset (in percentage) of the intervals to use for benchmarking. Format should be either XX-YY or XX,YY-ZZ. If this is left empty, all subsets/deciles are used |
| -b, --benchname | name of the benchmark, used for the output folder. This allows to perform multiple benchmarks |
| -c, --simname | name of the simulation data that is being used. Note that this should be the same as the name of the simulation data that was used for the simulation |
| -t, --tool | tool to benchmark. Currently, the following tools are supported: `tabix`, `bedtools`, `bedtools_sorted`, `bedtools_tabix` (deprecated), `bedops`, `bedmaps`, `giggle`, `granges`, `gia`, `bedtk`, `bedtk_sorted`, `igd`, `ailist`, `ucsc`, `awk`, `intervaltree` |

#### Benchmarked tools

Every command of the index and query steps is measured and summed (time) or maximised (memory); work that only prepares input or
converts the output into BED for the scoring is not measured.

| `--tool` | Index step (measured) | Query step (measured) | Notes |
| --- | --- | --- | --- |
| `tabix` | `sort`, `bgzip`, `tabix -C -p bed` (index size: `.gz` + `.csi`) | `tabix REF.bed.gz -R QUERY` | |
| `bedtools` | | `bedtools intersect -wa -a REF -b QUERY` | |
| `bedtools_sorted` | `sort` of the reference (no separate index) | `sort` of the query, `bedtools intersect -sorted -g GENOME -wa -a REF_SORTED -b QUERY_SORTED` | the sweep algorithm for sorted input; `-g` gives the chromosome order of the sorted data, so chromosomes present in only one file are handled; the genome file is written unmeasured |
| `bedtools_tabix` | as `tabix` | as `bedtools_sorted`, reading the bgzipped reference | deprecated, removed in 0.15.0: bedtools cannot use the tabix index for random access, so this measures `bedtools_sorted` plus an index cost |
| `bedops` | `sort` of the reference (no separate index) | `sort` of the query, `bedops --element-of 1 REF_SORTED QUERY_SORTED`; complex queries: `bedmap --echo-map --multidelim '\n' QUERY_SORTED REF_SORTED` | |
| `bedmaps` | `sort` of the reference (no separate index) | `sort` of the query, `bedmap --echo-map --multidelim '\n' QUERY_SORTED REF_SORTED` | |
| `giggle` | `giggle/scripts/sort_bed`, `giggle index -s` | `sort_bed` of the query, `giggle search -v` | |
| `granges` | | `granges filter --genome GENOME --left REF_SORTED --right QUERY` | reads the sorted reference; the `.tsv` copies and a genome file in natural chromosome order (granges 0.2.2 labels its query trees in that order, [#36](https://github.com/ylab-hi/segmeter/issues/36)) are prepared unmeasured |
| `gia` | | `gia intersect -a QUERY -b REF -t` | |
| `bedtk` | | `bedtk flt QUERY REF` | bedtk reports each reference interval once; the duplicates that complex queries expect are restored with an unmeasured `bedtools intersect` pass |
| `bedtk_sorted` | `sort` of the reference (no separate index) | `sort` of the query, `bedtk flt QUERY_SORTED REF_SORTED` | bedtk does not need sorted input, so this only adds the sorting cost |
| `igd` | `igd create` | `igd search -q QUERY -f` | the output is converted to BED unmeasured |
| `ailist` | | `ailist REF QUERY` | the overlap counts are expanded to one line per overlap unmeasured |
| `ucsc` | | `bedIntersect -aHitAny REF QUERY OUT` | |
| `awk` | | `tools/intersect_awk.py` | an awk script that hashes the intervals by chromosome, then scans that chromosome's intervals linearly |
| `intervaltree` | | `tools/intersect_intervaltree.py` | a Python script that builds one interval tree per chromosome (the `intervaltree` package) and queries it for every interval of the other file |

Not benchmarked: `gia intersect --sorted`. In gia 0.2.23 it numbers the chromosomes of each file by their order of appearance, so when one
file lacks a chromosome that the other has, every later chromosome is compared with the wrong one and overlaps are silently dropped
([gia#120](https://github.com/noamteyssier/gia/issues/120), [#53](https://github.com/ylab-hi/segmeter/issues/53)).

This generates a separate output folder for each benchmark tool in the folder `DATADIR/bench/benchname/` with a subfolder for each INTVLNUM.
In additional subfolders (`precision` and `stats`), the precision and statistics are stored.
The precision is stored in a file `DATADIR/bench/benchname/precision/<INTVLNUM>_<PERCENT>.txt` and the
statistics in `DATADIR/bench/benchname/stats/<INTVLNUM>_<PERCENT>.txt`.

With `--target` and `--query` (without `-r`), segmeter copies the target to `DATADIR/ref/target.bed` and prepares
`target_sorted.bed` and `target_chromlens.txt` (the chromosome lengths, needed by `granges`) next to it; this preparation
is not measured, like the sorted reference of the simulated data. The output folder `DATADIR/bench/benchname/TOOL/` then
contains `query_stats.txt` (time and memory of the query), `index_stats.txt` for index-based tools, `log.txt` with the
executed commands, and `result.bed` with the target intervals the tool reported as overlapping the query.

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

In addition, we provide a ready-to-use Docker container that has segmeter preconfigured. It can be found at [dockerhub](https://hub.docker.com/r/yanglabinfo/segmeter). We provide three different containers that can be used for the different tools.

| Container      | Tools      | Container tag |
| ------------- | ------------- | ------------- |
| giggle | giggle | segmeter:giggle-latest |
| others | ailist, bedops, bedtools, bedtk, igd, tabix, ucsc, awk, intervaltree | segmeter:others-latest |
| rust-tools | gia, granges | segmeter:rust-tools-latest |

The `latest` tags follow the most recent release; for reproducible runs use a version tag instead, e.g. `others-v0.13.2` (see [Published benchmark](#published-benchmark)).
The images are built for `linux/amd64` only (the arm64 images up to v0.14.0 ship an x86_64 BEDOPS). On Apple Silicon run them with
`--platform linux/amd64`, which Docker Desktop translates with Rosetta; such runs are fine to verify the setup, but their timings are not comparable to native ones.

This can be used with the following commands:
```
docker run -it -d -v /folder/on/host/:/folder/in/container/ yanglabinfo/segmeter:<container_tag> /bin/bash
docker exec <container_id> segmeter <args>
```

## Singularity

The Docker images can also be pulled with Singularity, e.g. `singularity pull docker://yanglabinfo/segmeter:others-v0.13.2`, which creates the `segmeter_others-v0.13.2.sif` file.
Consequently, this can be used with `singularity exec segmeter_others-v0.13.2.sif segmeter <args>`. Use the tag of the container that holds the tool you want to benchmark (see the table above).

## Published benchmark

segmeter is described in

> Schäfer RA, Yang R. A comprehensive benchmark of tools for efficient genomic interval querying. *Briefings in Bioinformatics*. 2025;26(4):bbaf379. [doi:10.1093/bib/bbaf379](https://doi.org/10.1093/bib/bbaf379)

The results in the article were produced with **segmeter v0.13.x** ([tags](https://github.com/ylab-hi/segmeter/tags)) and the following setup.
The tools were run in the three containers listed in the [Docker](#docker) section (`others`, `giggle`, `rust-tools`), as not all tools
build in the same environment. Use the images of the last patch release of v0.13.x, currently `yanglabinfo/segmeter:others-v0.13.2`,
`yanglabinfo/segmeter:giggle-v0.13.2` and `yanglabinfo/segmeter:rust-tools-v0.13.2`.

| Item | Setup |
| --- | --- |
| Simulated data | [doi:10.5281/zenodo.14880992](https://doi.org/10.5281/zenodo.14880992), created with `segmeter sim -o simdata -n 10,100,1K,10K,100K -c sim_001` |
| Simulation parameters | defaults: interval size 100-10000 bp (`-i`), gap size 100-5000 bp (`-g`), maximum chromosome length 1000000000 bp (`-m`) |
| Benchmark | `segmeter bench -o simdata -n 10,100,1K,10K,100K -c sim_001 -b bench_001 -t TOOL`, repeated three times (`bench_001`, `bench_002`, `bench_003`) with all subsets (`-s 10-100`, the default) |
| `--tool` choices in v0.13.x | `tabix`, `bedtools`, `bedtools_sorted`, `bedtools_tabix`, `bedops`, `bedmaps`, `giggle`, `granges`, `gia`, `bedtk`, `bedtk_sorted`, `igd`, `ailist`, `ucsc`, `awk`, `intervaltree` |
| `bedtools_sorted`, `bedtools_tabix`, `bedtk_sorted` in v0.13.x | run without `-sorted` on the simulator's sorted reference (`bedtools_tabix` with the unsorted query), the bgzip/tabix output of the index step unused: the v0.13.x index time of all three includes sort and bgzip (plus tabix for `bedtools_tabix`) and `index_size(MB)` is the bgzip (plus `.csi`) size. From 0.14.0 the variants read the index step's output, the bedtools variants pass `-sorted -g`, the `_sorted` variants index with the sort alone (index size reported as 0, like `bedops`), and `bedtools_tabix` is deprecated (its query step equals `bedtools_sorted`) ([#39](https://github.com/ylab-hi/segmeter/issues/39)) |
| `granges` in v0.13.x | reads the simulator's `<label>_chromlens.txt` as its genome file in the simulator's random draw order. granges 0.2.2 labels its query trees in natural chromosome order (1, 2, ..., 22, X, Y) but reads the genome file in file order, so every run whose genome file is not in natural order compared the left ranges of one chromosome with the queries of another: in the published `sim_001` data that is `1K`, `10K` and `100K` (`10` has one chromosome, `100` happened to be drawn in natural order), and their v0.13.x granges precision and recall are an artifact of the genome-file order, not of its overlap detection. From 0.14.0 granges reads an unmeasured copy of the genome file in natural order; the simulated data is unchanged ([#36](https://github.com/ylab-hi/segmeter/issues/36)) |
| Tool versions | bedtools 2.30.0, tabix (htslib) 1.16, BEDOPS 2.4.41, bedtk 0.0-r30, IGD 0.1.1, AIList 0.1.1, UCSC bedIntersect (kent source 482), intervaltree 3.1.0, GIGGLE 0.6.3, gia 0.2.23, granges 0.2.2 (as installed in the images above) |

To redo the benchmark, use the **last patch release of v0.13.x** (currently v0.13.2): it keeps these defaults, new options are off by default,
and the `bench` command takes the additional flag `-r`/`--simdata` for simulated data. Later minor releases (0.14.x and up) may change the tool versions and the simulated data.
