FROM ubuntu:18.04

RUN apt-get update && \
    apt-get -y install openjdk-8-jre-headless && \
    apt-get -y install curl wget && \
    apt-get -y install python3 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# we copy over the binaries under /home from the previous tag.
# and we include binaries from the go build
COPY fs/ /
