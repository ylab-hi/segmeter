# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# [Unreleased]
## Added
- Added tests for the measurement path: `tests/core/test_calls.py` checks that `tool_call` reports 0 MB when the `/usr/bin/time` output has no RSS line, and `tests/tools/test_query.py` runs every tool of `query_call` (with `index_call` for the index-based tools) on a random target/query pair and checks the reported reference intervals against a pure-Python oracle, skipping tools that are not installed; tests are split into `tests/core/` (tool-agnostic) and `tests/tools/` ([#12](https://github.com/ylab-hi/segmeter/issues/12), [#32](https://github.com/ylab-hi/segmeter/pull/32))
- Added optional `--max_span` (a multiple of 10) to `segmeter sim` to limit the number of intervals covered by a complex query, which otherwise makes the expected output grow quadratically with the intervals per chromosome; off by default, so the simulated data is unchanged ([#14](https://github.com/ylab-hi/segmeter/issues/14), [#23](https://github.com/ylab-hi/segmeter/pull/23))

## Fixed
- On macOS the memory of a tool call is converted from bytes: `/usr/bin/time -l` reports the resident set size in bytes, GNU time in kilobytes, and the fixed division by 1024 reported macOS memory about 1024 times too high (`true`: 1200 MB, now 1.17 MB). Linux and the containers, hence the published values, are unchanged ([#10](https://github.com/ylab-hi/segmeter/issues/10), [#55](https://github.com/ylab-hi/segmeter/pull/55))
- `tool_call` checks the exit code of the measured command (`/usr/bin/time` passes it through, 127 when the binary is missing) and raises with the code, the command and its stderr, which also lands in `log.txt` (closed on failure too): a tool that crashes, rejects its input (e.g. `granges filter` on BED3) or is not installed now aborts the benchmark instead of yielding a runtime, a memory value and an empty result. An empty query file (an empty complex bin, which occurs for small interval counts) runs no tool and reports 0 s / 0 MB, because `granges filter`, `gia intersect` and `gia sort` exit non-zero on an empty file and `tabix -R` dumps the whole reference; this changes the measured values of empty-bin rows for every tool (previously the tool's start-up on an empty file), their precision was already skipped ([#9](https://github.com/ylab-hi/segmeter/issues/9), [#54](https://github.com/ylab-hi/segmeter/pull/54))
- The `bedops` query folds the memory of its query sort into the reported query memory with `max`, as every other sorted variant does (its runtime was already counted). Changes the measured bedops query memory wherever the sort was the peak (the reported value was too low); runtime, precision and the other tools are unchanged. `tests/core/test_calls.py` runs every tool branch of `query_call` with a stubbed `tool_call` and checks that each step's runtime is added and its memory folded in ([#11](https://github.com/ylab-hi/segmeter/issues/11), [#51](https://github.com/ylab-hi/segmeter/pull/51))
- `query_call` no longer leaves temporary files behind: the intermediate files of a query (sorted query, raw output of bedtk/igd/ailist, granges' `.tsv` copies and genome file, giggle's sorted query) are written to one `tempfile.mkdtemp()` directory removed at the end of the call, and `BenchTool.query_intervals` unlinks the tool's output once it is scored (the real-data mode already moves it to `result.bed`); the tool commands are unchanged. `tests/tools/test_query.py` asserts after every tool that its temporary directory is empty ([#2](https://github.com/ylab-hi/segmeter/issues/2), [#50](https://github.com/ylab-hi/segmeter/pull/50))
- The `granges` query copies the reference and query into temporary files created with the `.tsv` suffix and closes them, instead of copying to `<name>.tsv` next to two empty temporary files that were never closed or removed (two files fewer per query in `$TMPDIR`; unlinking the remaining temporary files is #2), and `intersect_intervaltree.py` closes its output file when the query returns instead of at interpreter exit; the granges command and the tool output are unchanged ([#40](https://github.com/ylab-hi/segmeter/issues/40), [#49](https://github.com/ylab-hi/segmeter/pull/49))
- `granges` now reads a copy of the genome file in natural chromosome order (numbers, then X, Y, M), written unmeasured before the query: granges 0.2.2 labels its query trees in that order but reads the genome file in file order, so the simulator's random chromosome order made it compare the wrong chromosomes (a local 100K run scored precision 0.56 and recall 0.85 instead of 1.0). The simulated data is unchanged; the measured granges values with more than one chromosome change, so granges should be re-run. `tests/tools/test_query.py` writes the simulated chromosome lengths in reverse order to pin this ([#36](https://github.com/ylab-hi/segmeter/issues/36), [#46](https://github.com/ylab-hi/segmeter/pull/46))
- `bedtools_sorted`, `bedtools_tabix` and `bedtk_sorted` now query the output of their index step: the `_sorted` variants index with a sort only (no bgzip of an unread file) and read the sorted reference, `bedtools_tabix` reads the bgzipped, tabix-indexed reference, and both bedtools variants run `bedtools intersect -sorted -g` with the sorted query and a genome file in the order of the sorted data (`bedtools_tabix` used the unsorted query before); the `_sorted` variants report an index size of 0 like the other sort-only index steps. Changes the measured index and query values of these three tools; the README records the v0.13.x behavior ([#39](https://github.com/ylab-hi/segmeter/issues/39), [#42](https://github.com/ylab-hi/segmeter/pull/42))
- Fixed the real-data mode (`bench --target --query`): the target was queried against itself instead of the query file, and the sorted target and chromosome lengths that `bedtools_sorted`, `bedtools_tabix`, `bedmaps`, `bedtk_sorted` and `granges` read were never created, so they returned empty output; `BenchTool` now writes `target_sorted.bed` and `target_chromlens.txt` (in the chromosome order granges expects, numbers then X, Y, M), the simulation directory is no longer required in this mode, the tool's overlaps are kept as `bench/<benchname>/<tool>/result.bed`, and `bedops` reads that sorted target instead of sorting the reference inside the measured query (which its index step already measures) ([#19](https://github.com/ylab-hi/segmeter/issues/19), [#33](https://github.com/ylab-hi/segmeter/issues/33), [#37](https://github.com/ylab-hi/segmeter/pull/37))

## Refactored
- Removed the unreachable `gia_sorted` variant (never a `--tool` choice): gia 0.2.23's `intersect --sorted` numbers the chromosomes of each file by order of appearance and silently drops overlaps when one file lacks a chromosome of the other ([gia#120](https://github.com/noamteyssier/gia/issues/120)), so the variant is not benchmarked; no behavior change ([#53](https://github.com/ylab-hi/segmeter/issues/53), [#57](https://github.com/ylab-hi/segmeter/pull/57))
- Removed dead code found by a Python audit: unused simulator methods (`det_rightmost_start`, `update_leftgap`, `close_datafiles_complex`) and never-written file handles, the unused `num` argument of `index_call`/`query_call`, the never-read `reffiles["ref"]`, the unused `datatype` arguments of `get_query_group` and `sort_datafiles`, unread locals, imports and an unreachable guard, and the unused `--format`/`--stats` options of `intersect_intervaltree.py`; the scaffold branch of `select_chrom` is pinned by a test; no behavior change (seeded simulator output identical) ([#3](https://github.com/ylab-hi/segmeter/issues/3), [#38](https://github.com/ylab-hi/segmeter/pull/38))
- Cleanup from the Python audit: removed `SimBase` (a wrapper around `SimBED` behind the single-choice `--format`, which stays) and the `BenchTool.create_index`/`query_interval_file` wrappers, and the unreachable guard around `validate()`; one `utility.BASIC_QUERIES` table replaces three copies of the basic query types; `det_intvlnums` is one loop, `gapsize` is parsed once, the `igd` post-processing streams line by line; variables that shadowed builtins (`chr`, `dir`, `bin`, `sorted`) are renamed; `main.py` only runs `main()` when executed, so `tests/core/test_benchmark.py` can check `det_intvlnums`, `validate`, the negatives file and the recorded time and memory of `query_intervals`; no behavior change (seeded simulator output and bench results identical) ([#41](https://github.com/ylab-hi/segmeter/issues/41), [#48](https://github.com/ylab-hi/segmeter/pull/48))

## Changed
- README: new "Benchmarked tools" section with a table of every `--tool` (measured index step, measured query command, what is prepared or converted unmeasured), replacing the prose on the sorted variants, and a note on why `gia intersect --sorted` is not benchmarked: gia 0.2.23 numbers the chromosomes of each file by order of appearance and drops overlaps when one file lacks a chromosome of the other ([gia#120](https://github.com/noamteyssier/gia/issues/120)) ([#53](https://github.com/ylab-hi/segmeter/issues/53), [#56](https://github.com/ylab-hi/segmeter/pull/56))
- Deprecated `bedtools_tabix`: `bench -t bedtools_tabix` prints a warning on stderr, since the variant measures `bedtools_sorted` plus a tabix index that bedtools cannot use; the tool stays available until 0.15.0 ([#43](https://github.com/ylab-hi/segmeter/issues/43), [#44](https://github.com/ylab-hi/segmeter/pull/44))
- README: the "Published benchmark" table records the v0.13.x granges behavior (genome file in the simulator's random chromosome order, so the published granges precision and recall are an artifact of the genome-file order; corrected from 0.14.0) ([#36](https://github.com/ylab-hi/segmeter/issues/36), [#47](https://github.com/ylab-hi/segmeter/pull/47))
- README: the "Published benchmark" section moved to the end and refers to the last patch release of v0.13.x; Singularity instructions pull a versioned image instead of the stale unqualified `latest` tag; argument tables corrected (`--tool` lists all 16 tools, bench-only options moved to the bench table, `awk` in the container table) ([#31](https://github.com/ylab-hi/segmeter/pull/31))

# [0.13.2]
## Changed
- Pinned the tool versions in the containers to those of the published benchmark (bedtools 2.30.0, tabix 1.16, UCSC bedIntersect built from kent source 482, bedtk/IGD/AIList from their upstream repositories at the snapshot commits, giggle 0.6.3, gia 0.2.23 and granges 0.2.2 with Rust 1.87.0; `others` on `python:3.10-slim-bookworm`), build each image from the release tag it is named after, and replaced the deprecated `set-output` in the release workflows ([#17](https://github.com/ylab-hi/segmeter/issues/17), [#26](https://github.com/ylab-hi/segmeter/pull/26))
- Added the README section "Published benchmark" with the citation, the segmeter version (v0.13.x), the containers and tool versions, the Zenodo dataset, the commands and parameters of the published benchmark, and v0.13.2 as the release to redo it with ([#24](https://github.com/ylab-hi/segmeter/pull/24))
- Release workflows tag the giggle image by version and build it for arm64; the rust-tools image is tagged `rust-tools-latest` instead of `latest`; README documents `-r/--simdata`, `--query`, `--target` and the tools of each container image
- Removed unused `save_index_time` (`utility.py`) and a no-op memory comparison (`calls.py`); collapsed repeated max-memory blocks into `max()` calls
- Added `.gitignore` for `__pycache__`, `.pyc`, `.DS_Store`, `.Rhistory`

## Fixed
- Fixed the `segmeter` entry point in the containers, which failed with `from: command not found` because `main.py` had no shebang ([#25](https://github.com/ylab-hi/segmeter/issues/25), [#29](https://github.com/ylab-hi/segmeter/pull/29))
- Fixed quadratic precision scoring (tool results are looked up in a set) and count complex-query results without loading the output into memory; fixed the crash of `bench -r` since v0.13.1 (`querydirs` condition inverted when `--realdata` became `--simdata`) ([#13](https://github.com/ylab-hi/segmeter/issues/13), [#21](https://github.com/ylab-hi/segmeter/issues/21), [#22](https://github.com/ylab-hi/segmeter/pull/22))
- Fixed four correctness bugs in the measurement path: scaffold-name filter, bedops temp-file redirect, negative RSS sentinel, and unchecked bedtk dedup exit code ([#1](https://github.com/ylab-hi/segmeter/issues/1), [#7](https://github.com/ylab-hi/segmeter/pull/7))
- Fixed `file_linecounter` which always returned 1 instead of counting lines

# [0.13.1]
## Changed
- Renamed `--realdata` to `-r/--simdata`; added `--query` and `--target` for benchmarking arbitrary query/target files

## Fix
- Changed build container of rust-tools to new tag

# [0.13.0]
## Features
- Added support for intersection with awk
- Added support for intersection with interval tree (https://pypi.org/project/intervaltree/)

# [0.12.0]
## Features
- Added logfile to output in benchmark for each tool
- final version for manuscript

# [0.11.0]
## Features
- Added support for individual sim/bench runs (output defined by name)
- Renamed outputdir to datadir

# [0.10.0]
## Fix
- Benchmarking name can be specified (to allow multiple benchmarks)

# [0.9.2]
## Fix
- Changed base image in Github action for giggle to ubuntu:20.04

# [0.9.1]
## Fix
- Use internal sort function for gia (rather than unix sort)

# [0.9.0]
## Features
- Added support for GIA (unsorted and sorted)
- Added support for bedtk (unsorted and sorted)
- Added support for igd
- Added support for AIList
- Added support for gia (unsorted and sorted)
- Sort output of bench into label subfolders

# [0.8.0]
## Features
- Added support for granges (https://github.com/vsbuffalo/granges)

# [0.7.13]
## Fix
-  modify the Dockerfile to use the relative path from the root

# [0.7.12]
## Fix
- made some changes to the Dockerfile for giggle and its corresponding action

# [0.7.11]
## Fix
- context for giggle action within giggle container

# [0.7.10]
## Fix
- Adjust path in giggle index (it generates it in parent folder)

# [0.7.9]
## Fix
- another change in the context

# [0.7.8]
## Fix
- testing Dockerfile for dir content /segmeter

# [0.7.7]
## Fix
- added report in giggle action

# [0.7.6]
## Fix
- changed context of github action for giggle container - this ensures that only the segmeter path is added (rather than the whole repo path)

# [0.7.5]
## Fix
- added correct call to giggle

# [0.7.4]
## Fix
- Changed context for giggle docker action

# [0.7.3]
## Fix
- Fixed path in giggle docker action

# [0.7.2]
## Fix
- Fixed bug in Dockerfile for giggle

# [0.7.1]
## Fix
- Fixed bug in non-overlapping intervals
- Changed output format to be more readable

# [0.7.0]
## Features
- Separate subsets of the benchmark can be run using the `--subset` parameter

# [0.6.2]
## Fix
- Added missing Dockerfile for giggle-specific container

# [0.6.1]
## Fix
- It helps if the correct action is uploaded

# [0.6.0]
## Features
- added giggle to benchmark (in docker)

# [0.5.3]
## Fix
- Fixed wrong call in bedtools using random access with tabix

# [0.5.2]
## Fix
- fixed wrong call in bgzip to pipe the output

# [0.5.1]
## Fix
- f-string wrongly formatted caused file not found error

# [0.5.0]
## Feature
- added bedtools (with tabix and on sorted files) to benchmark
- code cleanup

# [0.4.0]

## Feature
- added bedops to benchmark

# [0.3.2]

## Fix
- reffiles subscriptable error fixed

# [0.3.1]

## Fix
- bedtools benchmark was missing the unsorted reference file

# [0.3.0]
- added bedtools to benchmark

# [0.2.9]

## Fix
- use parameter -R in tabix (for fairness) which specifies whole files

# [0.2.8]

## Fix
- fixed bug for single entries in intvlsize

# [0.2.7]

## Fix
- fixed bug with missing BenchTabix class

# [0.2.6]

## Fix
- fixed bug when only one number is provided (no multiple in comma separated list)

# [0.2.5]

- removed unnecessary print statements

# [0.2.4]

## Fix

- code refactoring for better modularization (finalize tabix benchmark)

# [0.2.3]

## Fix

- removed psutil

# [0.2.2]

## Fix

- removed wrong lib import

# [0.2.1]

## Fix

- added install instruction to time in Dockerfile

# [0.2.0]

## Feature

- minimum and maximum of the randomly generated gapsizes can be specified as paramter (--gapsize)
- minimum and maximum of the randomly generated interval sizes can be specified as paramter (--intvlsize)
- included memory measurements for tabix benchmark

# [0.1.2]

## Fix

- removed wrong import statements and logger

# [0.1.1]

## Fix

- added symlink in docker container and added missing files

# [0.1.0]

## Features

- Initial release
- Simulate 'simple' interval data
- Benchmark with tabix
