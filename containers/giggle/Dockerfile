FROM ubuntu:24.04
LABEL authors="Richard A. Schäfer"

WORKDIR /
RUN apt-get -y update
# giggle links against the system htslib since 2026 (libhts-dev, which brings the gnutls flavour of libcurl)
RUN apt-get install -y \
    gcc \
    make \
    autoconf \
    zlib1g-dev \
    libbz2-dev \
    libhts-dev \
    libssl-dev

RUN apt-get install -y \
    ruby \
    git \
    python3 \
    python3-pip \
    python3-dev \
    python3-setuptools \
    python3-wheel \
    time

# create symlink to python3
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN apt-get install -y \
    tabix

# pinned to upstream giggle (0.6.3 with the 2026 fixes). The published benchmark built the fork riasc/giggle at 1b7cb90,
# which removed giggle's check that the input and output directories share a parent; segmeter writes the index next to
# the input directory since #27, so upstream builds and runs unchanged
RUN git clone https://github.com/ryanlayer/giggle.git && cd /giggle && git checkout 215bf202d8e766c40f49faaf7f18dc19226a6d63
RUN cd /giggle && make HTS_INC=/usr/include HTS_LIB=/usr/lib/x86_64-linux-gnu

ENV GIGGLE_ROOT="/giggle" PATH="/giggle/bin/:$PATH"

WORKDIR /
# pin segmeter to the release that the image is built for (the workflow passes the tag)
ARG SEGMETER_REF=main
RUN git clone --branch ${SEGMETER_REF} --depth 1 https://github.com/ylab-hi/segmeter.git
# create symlink to segmeter/segmeter/main.py
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
