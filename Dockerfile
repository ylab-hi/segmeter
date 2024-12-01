############################################################
# Dockerfile to build the segmeter container
############################################################

FROM python:3.10-slim
LABEL authors="Richard A. Schäfer"

RUN apt-get -y update && apt-get install -y
RUN apt-get install -y \
    python3-dev \
    build-essential \
    zlib1g-dev \
    pkg-config \
    cmake \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    tree \
    curl \
    git \
    ca-certificates

# install rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# install tabix
RUN apt-get install -y tabix

# install bedtools
RUN apt-get install -y bedtools

#install granges
RUN cargo install granges

# load segmeter
ADD segmeter /segmeter
RUN chmod +x /segmeter/main.py
RUN ln -s /segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
