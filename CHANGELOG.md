# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
