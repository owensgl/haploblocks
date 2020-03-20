FROM ubuntu:14.04
LABEL maintainer="Jean-Sebastien Legare <jslegare@cs.ubc.ca>"

# This container will build the go sources found in the genomics
# repository. It will expect the sources repository to be
# mounted in /build/src/genomics, and will place the output in
# /build/bin
#
# See accompanied entrypoint.sh for details.


RUN apt-get update && \
    apt-get -y install software-properties-common make git && \
    add-apt-repository ppa:gophers/archive && \
    apt-get update && \
    apt-get -y install golang-1.9-go && \
    apt-get clean && \
    mkdir -p \
    	  /build/src /build/pkg /build/bin \
	  /usr/local/go/src /usr/local/go/pkg /usr/local/go/bin \
	  /go/src /go/pkg /go/bin \
	  /usr/local/bin

# Go build dependencies
ENV PATH=/usr/local/bin:/usr/lib/go-1.9/bin:/go/bin/:/usr/local/go/bin:/build/bin/:${PATH}
ENV GOPATH=/usr/local/go
RUN go get -u github.com/golang/lint/golint && \
    go get -u github.com/golang/dep/cmd/dep
ENV GOPATH=/go:/usr/local/go:/build

COPY entrypoint.sh dumb-init_1.2.1 /usr/local/bin/
ENTRYPOINT ["/usr/local/bin/dumb-init_1.2.1", "--", "/usr/local/bin/entrypoint.sh"]
# default parameters for entrypoint
CMD ["--build"]
