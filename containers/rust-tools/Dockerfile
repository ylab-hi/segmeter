############################################################
# Dockerfile to build the segmeter container
# This is only for rust-based tools: granges, gia
############################################################

FROM ubuntu:22.04
LABEL authors="Richard A. Schäfer"

RUN apt-get -y update && apt-get install -y
RUN apt-get install -y \
    build-essential \
    zlib1g-dev \
    # pkg-config \
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

# install rust
# RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --no-modify-path --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"

RUN cargo install granges # install granges
RUN cargo install gia # install gia

RUN git clone https://github.com/ylab-hi/segmeter.git
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
