############################################################
# Dockerfile to build the segmeter container
############################################################

FROM python:3.10-slim-trixie
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

# install tabix and bedtools
# pinned (Debian trixie); the published benchmark (v0.13.x) used tabix 1.16 and bedtools 2.30.0 of bookworm (#27).
# A trixie point release of either package removes the pinned version from the archive and fails this line: bump the pin
RUN apt-get install -y tabix=1.21+ds-1 bedtools=2.31.1+dfsg-2

# install bedops by link (the published benchmark used 2.4.41)
RUN curl -L -o  bedops.tar.bz2 https://github.com/bedops/bedops/releases/download/v2.4.42/bedops_linux_x86_64-v2.4.42.tar.bz2
RUN mkdir bedops
RUN tar -xvf bedops.tar.bz2 -C bedops
RUN cp bedops/bin/* /usr/local/bin/

# install intervaltree library (the published benchmark used 3.1.0)
RUN pip install intervaltree==3.2.1

# bedtk is pinned to release 1.2 (r34); IGD and AIList to the upstream commits of the published benchmark
# (identical to the snapshots in tools-src/ of 2025-01-25, unchanged upstream since); AIList is the 2019 version,
# later upstream versions changed the command line and the output format
RUN git clone https://github.com/lh3/bedtk.git && cd /bedtk && git checkout fa2cc15a8ab61d28bf6a4cfd50f36f883021a385
RUN cd /bedtk && make
RUN cp /bedtk/bedtk /usr/local/bin/

RUN git clone https://github.com/databio/IGD.git && cd /IGD && git checkout 4197c239380059e7e3d6eafce10011fc96047f5a
# igd create prints a progress dot every n_files/10 files, which divides by zero for fewer than 10 input files (segmeter
# indexes one file per label); the gcc 12 build of the published image compiled the modulo away, the gcc 14 build traps
# with SIGFPE, so the divisor is floored at 1 (a progress print, no effect on the index or the measured values), #27
RUN cd /IGD && sed -i "s|nf10 = n_files/10;|nf10 = n_files/10 > 0 ? n_files/10 : 1;|; s|nf10=n_files/10;|nf10 = n_files/10 > 0 ? n_files/10 : 1;|" src/igd_create.c && grep -c "n_files/10 : 1" src/igd_create.c
RUN cd /IGD && make
RUN cp /IGD/bin/igd /usr/local/bin/

RUN git clone https://github.com/databio/AIList.git && cd /AIList && git checkout d7fcddcd08668352e8e2836dae1c9bc622496a3d
RUN cd /AIList && make
RUN cp /AIList/bin/ailist /usr/local/bin/

# build UCSC bedIntersect from the archived kent source (version 502; the published benchmark used 482);
# UCSC publishes only unversioned binaries (and x86_64 only), the sources are archived per version. The v503 and
# v504 archives are not self-contained (the zlib-ng submodule sources and checkCompileFlags.sh are missing), #27
RUN apt-get install -y libpng-dev libssl-dev uuid-dev default-libmysqlclient-dev libbz2-dev libcurl4-openssl-dev liblzma-dev
RUN curl -L -o userApps.v502.src.tgz http://hgdownload.soe.ucsc.edu/admin/exe/userApps.archive/userApps.v502.src.tgz
RUN tar -xzf userApps.v502.src.tgz
RUN cd userApps/kent/src && MACHTYPE=local NOSQLTEST=1 make topLibs hgLib && MACHTYPE=local NOSQLTEST=1 make -C hg/bedIntersect compile
RUN cp userApps/kent/src/hg/bedIntersect/bedIntersect /usr/local/bin/

# load segmeter
# Build the checked-out source, including PR changes; release workflows check out the tag.
COPY segmeter /segmeter/segmeter
COPY containers/tool-versions.json /opt/segmeter-tool-versions.json
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
