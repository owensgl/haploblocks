SHELL := /bin/bash

# The name of the executables (default is current directory name)
BUILD_TARGETS := \
	align/align \
	reduce/reduce \
	vc/vc \
	cas/cas

PACKAGE = genomics
PACKAGE_INST = \
	$(PACKAGE)/align \
	$(PACKAGE)/reduce \
	$(PACKAGE)/vc \
	$(PACKAGE)/cas

.DEFAULT_GOAL: $(BUILD_TARGETS)

# These will be provided to the target
VERSION := $(shell cat VERSION)
BUILD := $(shell git rev-parse HEAD 2>/dev/null || cat ./BUILDNO)

# Use linker flags to provide version/build settings to the target
LDFLAGS=-ldflags "-X=main.Version=$(VERSION) -X=main.Build=$(BUILD)" $(XLDFLAGS)

# go source files, ignore vendor directory
SRC = $(shell find . -type f -name '*.go' -not -path "./vendor/*")

.PHONY: all build clean install fmt simplify check dep

all: build

$(BUILD_TARGETS): $(SRC)
	@echo "Building $(PACKAGE)/$(dir $@)"
	cd $(dir $@) && go build $(XBUILDFLAGS) $(LDFLAGS)

build: $(BUILD_TARGETS)
	@true

dep:
	: update the ./vendor directory with dependencies found
	dep ensure

clean:
	rm -f $(BUILD_TARGETS)

install:
	go install $(LDFLAGS) $(PACKAGE_INST)

#uninstall: clean
#	@rm -f $$(which ${TARGET})

fmt:
	gofmt -l -w $(SRC)

simplify:
	gofmt -s -l -w $(SRC)

check:
	@for d in $(dir $(BUILD_TARGETS)); do echo "checking '$$d'..."; test -z $(shell gofmt -l "$$d/main.go" | tee /dev/stderr) || echo "[WARN] Fix formatting issues with 'make fmt'"; done
	@for d in $$(go list ./... | grep -v /vendor/ | sed s@^_/@/@g); do golint "$${d}"; done
	@go tool vet ${SRC}

#run: install
#	@$(TARGET)
