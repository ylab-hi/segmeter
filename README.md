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




## Docker

In additon, we provide a ready-to-use Docker container that has segmeter preconfigured.

## Singularity

Segmeter is provided as docker container that can be used using `docker pull yanglabinfo/segmeter`.
