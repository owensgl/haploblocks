FROM ubuntu:18.04
LABEL maintainer="Jean-Sebastien Legare <jslegare@cs.ubc.ca>"

RUN apt-get update && \
    apt-get -y install software-properties-common make git && \
    add-apt-repository ppa:gophers/archive && \
    apt-get update && \
    apt-get -y install golang-1.9-go && \
    apt-get clean

RUN mkdir -p \
     	  /build/src /build/pkg /build/bin \
 	  /usr/local/go/src /usr/local/go/pkg /usr/local/go/bin \
 	  /go/src /go/pkg /go/bin \
 	  /usr/local/bin

# RUN apk update && \
#     apk add git && \
#     apk add bash && \
#     apk add make

# Go build dependencies
ENV PATH=/usr/local/bin:/usr/lib/go-1.9/bin:/go/bin/:/usr/local/go/bin:/build/bin/:${PATH}
ENV GOPATH=/usr/local/go
RUN go get -u github.com/golang/lint/golint && \
    go get -u github.com/golang/dep/cmd/dep


#COPY entrypoint.sh dumb-init_1.2.1 /usr/local/bin/
#ENTRYPOINT ["/usr/local/bin/dumb-init_1.2.1", "--", "/usr/local/bin/entrypoint.sh"]
# default parameters for entrypoint
#CMD ["--build"]
