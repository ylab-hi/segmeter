FROM ubuntu:20.04
LABEL authors="Richard A. Schäfer"

WORKDIR /
RUN apt-get -y update
RUN apt-get install -y \
    gcc \
    make \
    autoconf \
    zlib1g-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
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

# pinned to the commit of the published benchmark (giggle 0.6.3)
RUN git clone https://github.com/riasc/giggle.git && cd /giggle && git checkout 1b7cb90a65e89f4f2bfb8351b8191339acfc5616
RUN cd /giggle && make

ENV GIGGLE_ROOT="/giggle" PATH="/giggle/bin/:$PATH"

WORKDIR /
# pin segmeter to the release that the image is built for (the workflow passes the tag)
ARG SEGMETER_REF=main
RUN git clone --branch ${SEGMETER_REF} --depth 1 https://github.com/ylab-hi/segmeter.git
# create symlink to segmeter/segmeter/main.py
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
