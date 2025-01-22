FROM ubuntu:18.04
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
    git

RUN git clone https://github.com/riasc/giggle.git \
    && cd giggle \
    && make

ENV GIGGLE_ROOT="/giggle" PATH="/giggle/bin/:$PATH"

# now add segmeter
ADD ../segmeter /segmeter
RUN chmod +x /segmeter/main.py
RUN ln -s /segmeter/main.py /usr/local/bin/segmeter && chmod +x /usr/local/bin/segmeter
