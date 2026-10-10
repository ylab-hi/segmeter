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
segmeter sim -o DATADIR [-h] [-n INVLNUMS] [-m MAX_CHROMLEN] [-c SIMNAME] [-g GAPSIZE] [-i INTVLSIZE] [--max_span MAX_SPAN] [--seed SEED]

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
| --seed | seed of the random simulation. The same seed and parameters give the same data; a random seed is drawn by default. The seed used is written with the other parameters to `DATADIR/sim/simname/BED/<INTVLNUM>_parameters.txt` next to the data of each size, so any simulation can be repeated; the record is replaced with the data when a size is simulated again and kept when another size is added later. The sizes consume the random stream in order, so `-n 10K,100K --seed S` reproduces the `10K` files of `-n 10K --seed S`; to reproduce a size from a multi-size run, rerun the recorded `-n` list |

This will generates output files in the 'DATADIR/sim/simname/BED' folder. The files are in BED format (currently the only supported format) and can be used for benchmarking. In particular, the filers are located in the following folders:

```
DATADIR/sim/simname/BED/ref/ # reference intervals
DATADIR/sim/simname/BED/basic/ # basic queries
DATADIR/sim/simname/BED/complex/ # complex queries
```

In addition, segmeter generates for each specified `INTVLNUM`, a file with the length of each simulated chromosome (`DATADIR/sim/simname/BED/<INTVLNUM>_chromlens.txt`),
and the number of intervals per chromosome (`DATADIR/sim/simname/BED/<INTVLNUM>_chrnums.txt`), and the seed and the parameters of the run that produced it (`DATADIR/sim/simname/BED/<INTVLNUM>_parameters.txt`).

#### Reference

`DATADIR/sim/simname/BED/ref` contains the reference intervals. This is a BED4 file with the interval ID in the fourth column.
For each value in `INTVLNUMS`, there is a corresponding file with the intervals.

#### Basic queries

`DATADIR/sim/simname/BED/basic` contains the basic queries. In total, segmeter generates ten basic queries for each interval in the reference.
For each value in `INTVLNUMS`, there is a corresponding folder with the basic queries. In each folder, there is a  subfolders for the different
types of basic queries:
```
DATADIR/sim/simname/BED/basic/perfect/ # perfect overlaps (boundaries are identical with reference)
DATADIR/sim/simname/BED/basic/5p-partial/ # partial overlaps on 5' end
DATADIR/sim/simname/BED/basic/3p-partial/ # partial overlaps on 3' end
DATADIR/sim/simname/BED/basic/contained/ # overlap is contained within the reference
DATADIR/sim/simname/BED/basic/enclosed/ # overlap encloses the reference
DATADIR/sim/simname/BED/basic/perfect-gap/ # perfect overlap with a gap between reference intervals
DATADIR/sim/simname/BED/basic/left-adjacent-gap/ # adjacent to the 5'-end of the interval (no overlap)
DATADIR/sim/simname/BED/basic/right-adjacent-gap/ # adjacent to the 3'-end of the interval (no overlap)
DATADIR/sim/simname/BED/basic/mid-gap1/ # random overlap with a gap
DATADIR/sim/simname/BED/basic/mid-gap2/ # random overlap with a gap
```

In each of the subfolders, there is a BED4 file for each of the specified INTVLNUMS (e.g., `DATADIR/sim/simname/BED/basic/query/<query_type>/<INTVLNUM>.bed`).
In addition, the queries are subsample to 10-100% of the queries and stored in corresponding files (e.g., `DATADIR/sim/simname/BED/basic/query/<query_type>/INTVLNUM_<PERCENT>p.bed`).

##### Truth

In addition, the truth files are stored in `DATADIR/sim/simname/BED/basic/truth/<INTVLNUM>.bed`. These files contain the queries and their corresponding reference intervals. In each line
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

`DATADIR/sim/simname/BED/complex` contains the complex queries. Currently, this only includes `mult` queries which basically cover multiple reference intervals. For each chromosome, there is one query
for each number of covered intervals, from 2 up to the number of intervals on the chromosome or `--max_span`, whichever is smaller. According to the number of intervals that
are covered in a complex query, the queries are stored in deciles (e.g., `DATADIR/sim/simname/BED/complex/query/mult/<INTVLNUM>_<DECILE>bin.bed`). Again this contains the queries in BED4 format with and
identifier in the fourth column. Note that `mult_13` indicates that this query covers 13 reference intervals:
```
chr12	45837	160962	mult_13
chr12	146101	264093	mult_14
chr12	87719	214921	mult_15
chr12	122976	259591	mult_16
chr14	20381	105691	mult_13
```

##### Truth

The truth files for the complex queries are stored in `DATADIR/sim/simname/BED/complex/truth/<INTVLNUM>.bed`. This consists of the query interval and the corresponding number of intervals that are covered by the query.
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
| -t, --tool | tool to benchmark. Currently, the following tools are supported: `tabix`, `bedtools`, `bedtools_sorted`, `bedops`, `bedmaps`, `giggle`, `granges`, `gia`, `bedtk`, `bedtk_sorted`, `igd`, `ailist`, `ucsc`, `awk`, `intervaltree` |

#### Benchmarked tools

Every command of the index and query steps is measured and summed (time) or maximised (memory); work that only prepares input or
converts the output into BED for the scoring is not measured. The unmeasured steps are named in the Notes column: the genome file of
`bedtools_sorted`, the `.tsv` copies and the reordered genome file of `granges`, the shrunk query of `giggle`, the conversion of the `igd`
and `ailist` output into BED lines, and the `bedtools intersect` pass that restores the duplicates of `bedtk`, `granges` and `ucsc` for
the complex queries (see [Precision](#precision) for what that pass means for the complex score).

| `--tool` | Index step (measured) | Query step (measured) | Notes |
| --- | --- | --- | --- |
| `tabix` | `sort`, `bgzip`, `tabix -C -p bed` (index size: `.gz` + `.csi`) | `tabix REF.bed.gz -R QUERY` | |
| `bedtools` | | `bedtools intersect -wa -a REF -b QUERY` | |
| `bedtools_sorted` | `sort` of the reference (no separate index) | `sort` of the query, `bedtools intersect -sorted -g GENOME -wa -a REF_SORTED -b QUERY_SORTED` | the sweep algorithm for sorted input; `-g` gives the chromosome order of the sorted data, so chromosomes present in only one file are handled; the genome file is written unmeasured |
| `bedops` | `sort` of the reference (no separate index) | `sort` of the query, `bedops --element-of 1 REF_SORTED QUERY_SORTED`; complex queries: `bedmap --echo-map --multidelim '\n' QUERY_SORTED REF_SORTED` | |
| `bedmaps` | `sort` of the reference (no separate index) | `sort` of the query, `bedmap --echo-map --multidelim '\n' QUERY_SORTED REF_SORTED` | |
| `giggle` | `giggle/scripts/sort_bed`, `giggle index -s` | `sort_bed` of the query, `giggle search -v` | giggle treats the indexed and the query intervals as closed on both ends, so a reference that merely touches a query would be a hit; the query is shrunk to `[start+1, end-1]` unmeasured, which gives the half-open result for every query of 2 bp or more (a 1 bp query becomes the point `[start, start]`, which still hits a reference ending at `start`, since a closed interval cannot be empty). The simulated data is unaffected either way: a gap query ends one coordinate before the next interval, so it never touches a reference, and the published giggle precision has no false positives; the shrink matters for `--target --query` data ([#70](https://github.com/ylab-hi/segmeter/issues/70)) |
| `granges` | | `granges filter --genome GENOME --left REF_SORTED --right QUERY` | reads the sorted reference; the `.tsv` copies and a genome file in natural chromosome order (granges 0.2.2 labels its query trees in that order, [#36](https://github.com/ylab-hi/segmeter/issues/36)) are prepared unmeasured; `filter` reports each reference once, the duplicates that complex queries expect are restored with an unmeasured `bedtools intersect` pass, as for `bedtk` ([#69](https://github.com/ylab-hi/segmeter/issues/69)) |
| `gia` | | `gia intersect -a QUERY -b REF -t` | |
| `bedtk` | | `bedtk flt QUERY REF` | bedtk reports each reference interval once; the duplicates that complex queries expect are restored with an unmeasured `bedtools intersect` pass, so the complex score of such tools checks the set of references they found and takes the pairs from bedtools |
| `bedtk_sorted` | `sort` of the reference (no separate index) | `sort` of the query, `bedtk flt QUERY_SORTED REF_SORTED` | bedtk does not need sorted input, so this only adds the sorting cost |
| `igd` | `igd create` | `igd search -q QUERY -f` | the output is converted to BED unmeasured |
| `ailist` | | `ailist REF QUERY` | the overlap counts are expanded to one line per overlap unmeasured |
| `ucsc` | | `bedIntersect -aHitAny REF QUERY OUT` | `-aHitAny` reports each reference once; the duplicates are restored with the unmeasured `bedtools intersect` pass, as for `bedtk`. Without `-aHitAny` bedIntersect prints one line per pair, but with the coordinates of the intersection, not of the reference, and a BED3 reference cannot be mapped back ([#69](https://github.com/ylab-hi/segmeter/issues/69)) |
| `awk` | | `tools/intersect_awk.py` | an awk script that hashes the intervals by chromosome, then scans that chromosome's intervals linearly, printing a reference once per query that hits it |
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

intvlnum	bin	TP	FP	FN	Precision	Recall	F1	distance
1000	10bin	900	0	0	1.0	1.0	1.0	0
```

The upper part of the file contains the precision, recall, and F1 score for the basic queries and subset (e.g., 10% of the queries).
The lower part scores the complex queries per decile (e.g., 10bin, the queries that cover up to 10% of the reference intervals per
chromosome). The tools print the reference interval of a hit, not the query, so the (query, reference) pairs cannot be told apart in the
output; what can be scored is the set of reference intervals reported for the whole bin and how often each one is reported:

- **TP, FP, FN, precision, recall, F1** are computed on the set of reference intervals: expected are the intervals the queries of the
  bin cover (a complex query runs from the start of one reference interval to the end of another, so it covers exactly the intervals
  inside it, which segmeter finds in the sorted reference), reported are the distinct intervals in the tool output. There is no TN,
  since a complex bin has no negative queries.
- **distance** scores the pairs: a tool has to print a reference interval once per query that covers it, so the expected count of an
  interval is the number of queries of the bin covering it; the distance is the sum over all intervals of the absolute difference
  between the reported and the expected count. A missed interval covered by three queries adds 3, an interval reported once too often
  adds 1, an interval reported although no query covers it adds 1; a missing and an extra pair do not cancel (the v0.13.x and v0.14.x
  distance was the difference of the two line totals, where they did). A perfect tool scores 0.

Most tools print one line per pair, because they answer query by query (`bedtools`, `tabix`, `bedmap`, `giggle`, `gia`, `igd`,
`intervaltree`; `ailist` prints a count per reference that segmeter expands into repeated lines). `bedtk flt`, `granges filter` and
`bedIntersect -aHitAny` answer the other way round and print each reference that is hit by any query once, so their output carries no
pairs (none of the three tools has a mode that prints the reference once per query: `bedtk isec` and plain `bedIntersect` print
intersections, `granges map` one aggregated line per reference). For them segmeter runs an unmeasured `bedtools intersect` pass over
their output and the complex query file: `-wa` prints each reported reference once per query that overlaps it, and `-v` appends the
reported references that no query overlaps, so the set scores see every reference the tool reported, and only the multiplicity of the
true hits comes from bedtools. For these three tools the distance therefore cannot show a pairing error of the tool itself. The pass
runs on the complex query files only; the basic queries are scored on the raw output of every tool. Checking the complex pairs
themselves would need one call per query, which is not feasible at benchmark scale.

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

### Tool versions

Every tool is pinned in its Dockerfile, so an image tag reproduces one set of versions. The current pins (images from 0.15.0) and the
versions of the published benchmark (images up to v0.14.x, see [Published benchmark](#published-benchmark)):

| Tool | Container | Current pin | Published benchmark | Pinned in |
| --- | --- | --- | --- | --- |
| bedtools | others | 2.31.1 (Debian trixie `2.31.1+dfsg-2`) | 2.30.0 (Debian bookworm) | `containers/others/Dockerfile`, apt |
| tabix (htslib) | others | 1.21 (Debian trixie `1.21+ds-1`) | 1.16 (Debian bookworm) | `containers/others/Dockerfile`, apt |
| BEDOPS (`bedops`, `bedmap`) | others | 2.4.42 | 2.4.41 | `containers/others/Dockerfile`, release tarball |
| bedtk | others | 1.2 (r34, commit `fa2cc15`) | 0.0-r30 (commit `da1fb73`, since tagged v1.0 upstream) | `containers/others/Dockerfile`, git commit |
| IGD | others | 0.1.1 (commit `4197c23`, 2021; built with a one-line patch that floors the divisor of a progress print at 1, since `igd create` divides by zero for fewer than 10 input files and gcc 14 no longer compiles that away) | same, unpatched (gcc 12) | `containers/others/Dockerfile`, git commit |
| AIList | others | 0.1.1 (commit `d7fcddc`, 2019; later versions changed the command line) | same | `containers/others/Dockerfile`, git commit |
| UCSC bedIntersect | others | kent source 502 | kent source 482 | `containers/others/Dockerfile`, archived source |
| intervaltree | others | 3.2.1 | 3.1.0 | `containers/others/Dockerfile`, pip |
| GIGGLE | giggle | 0.6.3, upstream commit `215bf20` (2026-03), Ubuntu 24.04, system htslib 1.19 | 0.6.3, fork `riasc/giggle` at `1b7cb90`, Ubuntu 20.04 | `containers/giggle/Dockerfile`, git commit |
| gia | rust-tools | 0.2.23, Rust 1.99.0 | 0.2.23, Rust 1.87.0 | `containers/rust-tools/Dockerfile`, cargo |
| granges | rust-tools | 0.2.2, Rust 1.99.0 | 0.2.2, Rust 1.87.0 | `containers/rust-tools/Dockerfile`, cargo |
| Python (harness) | all | 3.10 (others), 3.12 (giggle, Ubuntu 24.04), 3.10 (rust-tools, Ubuntu 22.04) | 3.10 (others), 3.8 (giggle, Ubuntu 20.04), 3.10 (rust-tools) | base images |

Measured values depend on the versions, so compare results only between runs of the same image tag. The pins are checked against
upstream before a release ([#27](https://github.com/ylab-hi/segmeter/issues/27)).

## Continuous benchmarking

The smoke workflow runs on pull requests and pushes to `main` in the three tool containers, with 10K and 100K intervals, seed 1729,
`--max_span 100`, and subset/bin 100. It checks every tool's basic and complex precision and recall, and requires complex
distance zero. Empty complex bins (in small full-run datasets) must have no false positives or missed hits.
Query timings are divided by bedtools' timing of the same query case on the same runner. PRs run the exact base commit
and the merge result sequentially using the same tool images; differing dataset hashes invalidate the comparison.
The comparison fails for a normalized slowdown above 25%. Smoke uses one sample per case, so inspect raw timings
and rerun a noisy failure before attributing it to a code change. Tool-image changes are evaluated on the full dashboard;
the PR comparison deliberately holds the tool environment constant.
Both tiers also record absolute simulation time and end-to-end time per tool (including scoring and container startup)
to detect regressions in segmeter outside the measured tool commands. The PR comparison checks these on the same runner too.

The full benchmark runs as a Slurm batch job using Singularity (or Apptainer), independently of GitHub Actions.
It covers 10, 100, 1K, 10K, 100K and 1M, all subsets/bins, with three independent index/query runs.
The exported absolute times are medians. Raw statistics, precision, SIF hashes, tool versions, dataset hashes,
a source snapshot, CPU information and Slurm allocation details are saved in the output directory.
Submit weekly or after a tool version changes; no GitHub runner is required on the cluster.

Both tiers export `customSmallerIsBetter` JSON for `benchmark-action/github-action-benchmark`. Main smoke runs and
full runs publish separate histories to `gh-pages` at `dev/bench/smoke`, `dev/bench/full-simulated` and `dev/bench/full-zenodo`. Enable GitHub Pages
for that branch to serve the dashboards. Use the same cluster node type and filesystem for comparable absolute timings. The full GitHub workflow only publishes
uploaded results; it does not submit cluster jobs.

To run locally (Docker required):

```sh
python3 scripts/continuous_benchmark.py --build --output /tmp/segmeter-smoke
python3 scripts/continuous_benchmark.py --build --tier full --output /tmp/segmeter-full
```

The output directory must be new. `--sizes 1K` permits a short integration check. Containers copy the local `segmeter/`
source at build time; the runner snapshots and mounts the selected source read-only at runtime, which lets it test both sides of a PR.
Update `containers/tool-versions.json` with any Dockerfile tool pin change; the manifest is embedded in each image and
read back into results. Timing comparisons across tool versions or different CPUs need separate interpretation.

### Full benchmark on a Slurm cluster

Run the launcher on your cluster's login/transfer node. It loads the module, pulls the three SIF containers,
optionally downloads and verifies the published Zenodo data, freezes the source and submits the benchmark to Slurm:

```sh
bash scripts/cluster/run-benchmark.sh \
  --module singularity/4.1 \
  --image-tag v0.14.1 \
  --work-dir /shared/project/segmeter-benchmarks \
  --dataset zenodo \
  --account YOUR_ACCOUNT --partition YOUR_PARTITION
```

Replace the module name, account, partition and shared path with your site's values, and choose a published image tag.
Use `--dataset simulated` (the default) for seeded data with the complex-span cap, or `--dataset zenodo` for the
published archive ([Zenodo record 14880992](https://zenodo.org/records/14880992)). The fixed archive is 216 MB compressed;
the launcher verifies its published MD5 checksum before extracting it. Zenodo contains sizes through 100K, so that run
omits 1M and records the archive provenance instead of claiming a seed or span cap. Published and seeded runs have
separate dashboard histories (`dev/bench/full-zenodo` and `dev/bench/full-simulated`).

For a setup check, append `--sizes 1K`. For resource overrides, add `--mem 64G`, `--time 2-00:00:00` or `--exclusive`.
Use `--module none` if the runtime is already on PATH; `--runtime apptainer --module YOUR_APPTAINER_MODULE` selects Apptainer.
Module loading is repeated inside the batch job. The launcher needs Python, Git, `flock`, network access and `sbatch` on
the login node. Compute nodes need Python and the runtime, but no network access or GitHub credentials.

Images and data are cached under the work directory. Each invocation prints the job ID, a unique results directory
and Slurm log path. Source is frozen before submission, so checkout edits while the job is queued do not change it.
Once the job finishes, `RESULTS/benchmark/publish/` contains the three JSON files to upload as described below.
Container pulls create the SIF files; `singularity exec` then runs them directly, without a separate persistent container.
Older release images are supported: when their embedded version manifest is absent, results record installed tool binary
SHA-256 fingerprints and package versions. They never reuse the current Dockerfile's pins as a claim about an old image.

The following commands are optional lower-level preparation and submission steps if you prefer to manage images yourself.

Use an x86_64 Linux node, Python 3.8+ and Git on the host, and SingularityCE 3.6+ (or Apptainer).
Load your site's Python and Singularity modules before submission; module names, account and partition are site-specific.
The job inherits that environment, then runs containers with a clean environment and explicit source/data binds
([Singularity exec options](https://docs.sylabs.io/guides/4.6/user-guide/cli/singularity_exec.html)).

First commit and push the benchmark code, then use a separate checkout on shared cluster storage. Keep that checkout
unchanged until the job starts; the job snapshots the Python source once running. The output and image directories
must be accessible from compute nodes. Local scratch can be used if you arrange to copy results back before allocation ends.

On a login/transfer node with internet access, pull a versioned release:

```sh
# Replace vX.Y.Z with the published release to benchmark; do not use latest.
bash scripts/cluster/pull-images.sh vX.Y.Z /shared/project/segmeter-images-vX.Y.Z
```

The new image directory's parent must exist. New images embed `/opt/segmeter-tool-versions.json`; older images use runtime binary fingerprints instead.
To benchmark the latest Dockerfile changes before release, build on a Docker machine and transfer Docker archives to the cluster:

```sh
# On the Docker machine, from this checkout; repeat for giggle and rust-tools.
docker build --platform linux/amd64 -t segmeter-bench:others -f containers/others/Dockerfile .
docker save --output others.tar segmeter-bench:others
# Transfer others.tar to the cluster, then on the login/transfer node:
singularity build /shared/project/segmeter-images/others.sif docker-archive://others.tar
```

Repeat the conversion for `giggle.sif` and `rust-tools.sif`. All three SIF files must be in the same image directory.
Pulling/conversion happens before submission; the compute job needs no internet access or GitHub credentials.
For Apptainer, use `SEGMETER_RUNTIME=apptainer` for both the pull script and `sbatch` invocation.

From the segmeter checkout, submit a small check first, then the full benchmark with a new output directory:

```sh
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/cluster/full-benchmark.sbatch \
  /shared/project/segmeter-images /shared/project/segmeter-check --sizes 1K

sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/cluster/full-benchmark.sbatch \
  /shared/project/segmeter-images /shared/project/segmeter-full-2026-10-09
```

The script requests one node, one task/CPU, 32 GB and 24 hours; these are starting resource requests, not measured
requirements. Override `--mem`, `--time`, and any site-specific resource options with `sbatch`. All tools run sequentially
within the allocation. Slurm writes `segmeter-full-JOBID.log` in the submission directory; monitor with `squeue`/`sacct`.
The 1M awk cases can take substantial time. For exclusive-node measurements, add `--exclusive` if your site permits it.
See [sbatch options](https://slurm.schedmd.com/sbatch.html) for resource overrides.

A successful full run creates three small files in `OUTPUT/publish/`; that directory is created only after every
correctness check passes. Keep the whole output directory and Slurm log on cluster storage for later inspection.
An abbreviated `--sizes` check is useful for setup validation but the publishing workflow rejects it.

### Publish cluster results

Copy `OUTPUT/publish/` to a machine with GitHub CLI access. The source commit must have been pushed to this repository.
Create a draft release to hold the results, upload the three JSON files, and dispatch the publishing workflow:

```sh
# Run from a segmeter checkout with gh authenticated; choose a unique results tag.
gh release create benchmark-results-2026-10-09 --draft --target MEASURED_COMMIT \
  --title "Cluster benchmark 2026-10-09" --notes "Slurm/Singularity full benchmark results"
gh release upload benchmark-results-2026-10-09 publish/segmeter-full-*.json
gh workflow run benchmark-full.yml -f results_release=benchmark-results-2026-10-09
```

`MEASURED_COMMIT` is `source_commit` in `segmeter-full-versions.json`. The workflow must already be on the default branch.
The results release can remain a draft. The workflow validates the files, retains them as an artifact, and publishes the
full dashboard to `gh-pages` under `dev/bench/full-simulated` or `dev/bench/full-zenodo`, attributed to the measured commit. Publishing requires repository write access and
GitHub Pages enabled for that branch. Cluster jobs and uploads are manual; there is no automatic weekly cluster submission.

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
| `bedtools_sorted`, `bedtools_tabix`, `bedtk_sorted` in v0.13.x | run without `-sorted` on the simulator's sorted reference (`bedtools_tabix` with the unsorted query), the bgzip/tabix output of the index step unused: the v0.13.x index time of all three includes sort and bgzip (plus tabix for `bedtools_tabix`) and `index_size(MB)` is the bgzip (plus `.csi`) size. From 0.14.0 the variants read the index step's output, the bedtools variants pass `-sorted -g`, the `_sorted` variants index with the sort alone (index size reported as 0, like `bedops`), and `bedtools_tabix` is deprecated (its query step equals `bedtools_sorted`), later removed ([#39](https://github.com/ylab-hi/segmeter/issues/39), [#45](https://github.com/ylab-hi/segmeter/issues/45)) |
| `granges` in v0.13.x | reads the simulator's `<label>_chromlens.txt` as its genome file in the simulator's random draw order. granges 0.2.2 labels its query trees in natural chromosome order (1, 2, ..., 22, X, Y) but reads the genome file in file order, so every run whose genome file is not in natural order compared the left ranges of one chromosome with the queries of another: in the published `sim_001` data that is `1K`, `10K` and `100K` (`10` has one chromosome, `100` happened to be drawn in natural order), and their v0.13.x granges precision and recall are an artifact of the genome-file order, not of its overlap detection. From 0.14.0 granges reads an unmeasured copy of the genome file in natural order; the simulated data is unchanged ([#36](https://github.com/ylab-hi/segmeter/issues/36)) |
| `awk`, `ucsc`, `granges` complex score in v0.13.x and v0.14.x | these tools report each reference once however many complex queries hit it (`bedIntersect -aHitAny`, `granges filter`, and the awk script stopped at the first hit), while the complex score counts one output line per (query, reference) pair; their complex distance is therefore the number of missing duplicates, not missed overlaps (the basic scores are unaffected, a basic query hits one reference). `bedtk` already had its duplicates restored by an unmeasured `bedtools intersect` pass, which ran on every query file, so a reference bedtk had reported for a basic gap query without overlapping it would have been dropped before the basic score (bedtk reported none). From 0.15.0 `granges` and `ucsc` get the same pass, it runs on the complex query files only, and the awk script prints a reference once per query ([#69](https://github.com/ylab-hi/segmeter/issues/69)) |
| Tool versions | bedtools 2.30.0, tabix (htslib) 1.16, BEDOPS 2.4.41, bedtk 0.0-r30, IGD 0.1.1, AIList 0.1.1, UCSC bedIntersect (kent source 482), intervaltree 3.1.0, GIGGLE 0.6.3, gia 0.2.23, granges 0.2.2 (as installed in the images above) |
| Tool versions from 0.15.0 | the images are pinned to the upstream versions of 2026-10, see the [tool version table](#tool-versions) in the Docker section. The measured values of every tool change with the versions, so results from 0.15.0 images are not comparable with the published ones ([#27](https://github.com/ylab-hi/segmeter/issues/27)) |

To redo the benchmark, use the **last patch release of v0.13.x** (currently v0.13.2): it keeps these defaults, new options are off by default,
and the `bench` command takes the additional flag `-r`/`--simdata` for simulated data. Later minor releases (0.14.x and up) may change the tool versions and the simulated data.
