FROM ubuntu:20.04 as builder
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy full-upgrade && \
    apt-get -qqy install --no-install-recommends \
    build-essential gcc  make zlib1g-dev \
    git ca-certificates
RUN git clone --depth 1 --branch v0.6 https://github.com/jsh58/Genrich.git
WORKDIR Genrich
RUN make

FROM ubuntu:20.04
COPY --from=builder /Genrich/Genrich /usr/local/bin/
WORKDIR /mnt

COPY run_tests.sh /usr/local/bin/
RUN  chmod 755 /usr/local/bin/run_tests.sh

ENTRYPOINT ["Genrich"]
CMD []


