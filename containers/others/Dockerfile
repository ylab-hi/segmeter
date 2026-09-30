############################################################
# Dockerfile to build the segmeter container
############################################################

FROM python:3.10-slim-bookworm
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
    ca-certificates \
    libc6-dev

# install tabix
# pinned to the versions of the published benchmark (Debian bookworm)
RUN apt-get install -y tabix=1.16+ds-3 bedtools=2.30.0+dfsg-3

# install bedops by link
RUN curl -L -o  bedops.tar.bz2 https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
RUN mkdir bedops
RUN tar -xvf bedops.tar.bz2 -C bedops
RUN cp bedops/bin/* /usr/local/bin/

# install intervaltree library
RUN pip install intervaltree==3.1.0

# install bedtk
RUN curl -L -o bedtk.tar.gz https://github.com/ylab-hi/segmeter/raw/refs/heads/main/tools-src/bedtk-v2025-01-25.tar.gz
RUN tar -xvf bedtk.tar.gz
RUN cd bedtk && make
RUN cp /bedtk/bedtk /usr/local/bin/

# install IGD
RUN curl -L -o IGD.tar.gz https://github.com/ylab-hi/segmeter/raw/refs/heads/main/tools-src/IGD-v2025-01-25.tar.gz
RUN tar -xvf IGD.tar.gz
RUN cd IGD && make
RUN cp /IGD/bin/igd /usr/local/bin/

# install AIList
RUN curl -L -o AIList.tar.gz https://github.com/ylab-hi/segmeter/raw/refs/heads/main/tools-src/AIList-v2025-01-25.tar.gz
RUN tar -xvf AIList.tar.gz
RUN cd AIList && make
RUN cp /AIList/bin/ailist /usr/local/bin/

# load UCSC tools (bedIntersect) - UCSC offers no versioned download; the published benchmark used kent source version 482
RUN curl -LO http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedIntersect
RUN chmod +x bedIntersect
RUN mv bedIntersect /usr/local/bin/

# load segmeter
# pin segmeter to the release that the image is built for (the workflow passes the tag)
ARG SEGMETER_REF=main
RUN git clone --branch ${SEGMETER_REF} --depth 1 https://github.com/ylab-hi/segmeter.git
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
