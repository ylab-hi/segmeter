# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
