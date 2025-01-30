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

RUN git clone https://github.com/riasc/giggle.git
RUN cd giggle
RUN cd lib/htslib
RUN autoreconf --install
RUN cd ../.. && make

ENV GIGGLE_ROOT="/giggle" PATH="/giggle/bin/:$PATH"

WORKDIR /
RUN git clone https://github.com/ylab-hi/segmeter.git
# create symlink to segmeter/segmeter/main.py
RUN ln -s /segmeter/segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
