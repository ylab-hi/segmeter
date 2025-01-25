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
    libssl-dev \
    ruby \
    git \
    python3 \
    python3-pip \
    python3-dev \
    python3-setuptools \
    python3-wheel \
    tabix \
    time

RUN git clone https://github.com/riasc/giggle.git \
    && cd giggle \
    ls -la \
    && make

ENV GIGGLE_ROOT="/giggle" PATH="/giggle/bin/:$PATH"

# now add segmeter
# ADD segmeter /segmeter
# RUN ls -l /segmeter
# RUN chmod +x /segmeter/main.py
# RUN ln -s /segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
