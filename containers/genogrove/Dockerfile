############################################################
# Dockerfile to build the segmeter container for genogrove, pygenogrove and polars-bio.
# genogrove needs GCC 13+ (C++20) and polars-bio needs Python 3.11+, which Ubuntu 24.04 provides (Python 3.12).
############################################################

FROM ubuntu:24.04
LABEL authors="Richard A. Schäfer"

RUN apt-get -y update && apt-get install -y \
    build-essential \
    cmake \
    pkg-config \
    git \
    ca-certificates \
    zlib1g-dev \
    libhts-dev=1.19+ds-1.1build3 \
    python3 \
    python3-pip \
    python3-venv \
    time

# genogrove CLI, pinned to a release tag; the CLI target is built in-tree (cxxopts is fetched by CPM at configure time)
ARG GENOGROVE_VERSION=v0.26.1
RUN git clone --branch ${GENOGROVE_VERSION} --depth 1 https://github.com/genogrove/genogrove.git /genogrove-src && \
    cmake -S /genogrove-src -B /genogrove-src/build -DCMAKE_BUILD_TYPE=Release -DBUILD_CLI=ON && \
    cmake --build /genogrove-src/build --target genogrove-cli -j"$(nproc)" && \
    cp /genogrove-src/build/bin/genogrove /usr/local/bin/genogrove && \
    rm -rf /genogrove-src

# pygenogrove (manylinux wheel) and polars-bio (manylinux abi3 wheel), pinned
RUN pip install --break-system-packages --no-cache-dir pygenogrove==0.9.0 polars-bio==0.36.2

WORKDIR /
# Build the checked-out source, including PR changes; release workflows check out the tag.
COPY segmeter /segmeter/segmeter
COPY containers/tool-versions.json /opt/segmeter-tool-versions.json
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
