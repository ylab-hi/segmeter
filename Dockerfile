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
    time \
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
# install gia
RUN cargo install gia

# install bedops by link
RUN curl -L -o  bedops.tar.bz2 https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
RUN mkdir bedops
RUN tar -xvf bedops.tar.bz2 -C bedops
RUN cp bedops/bin/* /usr/local/bin/

# install bedtk
RUN curl -L -o bedtk.tar.gz https://github.com/riasc/bedtk/releases/download/2025-01-27/bedtk-v2025-01-27.tar.gz
RUN tar -xvf bedtk.tar.gz
RUN cd bedtk && make
RUN cp /bedtk/bedtk /usr/local/bin/

# install IGD
RUN curl -L -o IGD.tar.gz https://github.com/riasc/IGD/releases/download/2025-01-27/IGD-2025-01-27.tar.gz
RUN tar -xvf IGD.tar.gz
RUN cd IGD && make
RUN cp /IGD/bin/igd /usr/local/bin/

# load segmeter
ADD segmeter /segmeter
RUN chmod +x /segmeter/main.py
RUN ln -s /segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
