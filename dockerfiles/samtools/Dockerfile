FROM ubuntu:20.04 as builder
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy full-upgrade && \
    apt-get -qqy install --no-install-recommends \
    git ca-certificates  build-essential gcc  make zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses-dev

RUN git clone --depth 1 https://github.com/samtools/samtools.git 
RUN git clone --depth 1 https://github.com/samtools/htslib.git 
WORKDIR /htslib
RUN git submodule update --init --recursive
RUN make && make test && make install
WORKDIR /samtools
RUN make && make test && make install

FROM ubuntu:20.04
COPY --from=builder /samtools/samtools /usr/local/bin/
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy full-upgrade && \
    apt-get -qqy install --no-install-recommends \
    libcurl4
WORKDIR /mnt
ENTRYPOINT ["samtools"]
CMD []
