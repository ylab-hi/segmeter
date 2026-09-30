############################################################
# Dockerfile to build the segmeter container
# This is only for rust-based tools: granges, gia
############################################################

FROM ubuntu:22.04
LABEL authors="Richard A. Schäfer"

RUN apt-get -y update && apt-get install -y
RUN apt-get install -y \
    build-essential
RUN apt-get install -y \
    zlib1g-dev
RUN apt-get install -y \
    pkg-config
RUN apt-get install -y \
    cmake

RUN apt-get install -y \
    python3-dev \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    tree \
    time \
    curl \
    git \
    ca-certificates \
    libc6-dev

# install rust - toolchain pinned to the one of the published benchmark image
# RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --no-modify-path --default-toolchain 1.87.0
ENV PATH="/root/.cargo/bin:${PATH}"

# pinned to the versions of the published benchmark
RUN cargo install granges --version 0.2.2 --locked
RUN cargo install gia --version 0.2.23 --locked

# pin segmeter to the release that the image is built for (the workflow passes the tag)
ARG SEGMETER_REF=main
RUN git clone --branch ${SEGMETER_REF} --depth 1 https://github.com/ylab-hi/segmeter.git
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
