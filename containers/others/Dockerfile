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

# bedtk, IGD and AIList are pinned to the upstream commits of the published benchmark
# (identical to the snapshots in tools-src/ of 2025-01-25); AIList is the 2019 version,
# later upstream versions changed the command line and the output format
RUN git clone https://github.com/lh3/bedtk.git && cd /bedtk && git checkout da1fb736ec6427d9021d9f586adbdd56aaf72330
RUN cd /bedtk && make
RUN cp /bedtk/bedtk /usr/local/bin/

RUN git clone https://github.com/databio/IGD.git && cd /IGD && git checkout 4197c239380059e7e3d6eafce10011fc96047f5a
RUN cd /IGD && make
RUN cp /IGD/bin/igd /usr/local/bin/

RUN git clone https://github.com/databio/AIList.git && cd /AIList && git checkout d7fcddcd08668352e8e2836dae1c9bc622496a3d
RUN cd /AIList && make
RUN cp /AIList/bin/ailist /usr/local/bin/

# build UCSC bedIntersect from the kent source of the published benchmark (version 482);
# UCSC publishes only unversioned binaries (and x86_64 only), the sources are archived per version
RUN apt-get install -y libpng-dev libssl-dev uuid-dev default-libmysqlclient-dev libbz2-dev
RUN curl -L -o userApps.v482.src.tgz http://hgdownload.soe.ucsc.edu/admin/exe/userApps.archive/userApps.v482.src.tgz
RUN tar -xzf userApps.v482.src.tgz
RUN cd userApps/kent/src && MACHTYPE=local NOSQLTEST=1 make topLibs hgLib && MACHTYPE=local NOSQLTEST=1 make -C hg/bedIntersect compile
RUN cp userApps/kent/src/hg/bedIntersect/bedIntersect /usr/local/bin/

# load segmeter
# pin segmeter to the release that the image is built for (the workflow passes the tag)
ARG SEGMETER_REF=main
RUN git clone --branch ${SEGMETER_REF} --depth 1 https://github.com/ylab-hi/segmeter.git
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
