#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# Build longranger.
#
TOP_DIR=$(shell pwd)
LIBPY=lib/python
LARIAT_DIR=lib/lariat
TADA_DIR=lib/tada
GOBINS=loupe/goloupe reads/bucket_fastq_by_bc reads/sort_fastq_by_bc cnv/molecular_count
VERSION=$(shell git describe --tags --always --dirty)

# Locations for python includes and libraries.
PYTHON_ROOT?=$(dir $(shell which python2.7))..
PYTHON_LIB_DIR?=$(dir $(lastword $(wildcard $(PYTHON_ROOT)/lib*/libpython2.7.so*)))
PYTHON_INC_DIR?=$(PYTHON_ROOT)/include/python2.7

export CFLAGS=-DLZ4_DISABLE_DEPRECATE_WARNINGS
export BUILDPATH=$(shell pwd)/lib/go
export GOPATH=$(shell pwd)/lib/go:$(shell pwd)/tenkit/lib/go

GO_VERSION=$(strip $(shell go version | sed 's/.*go\([0-9]*\.[0-9]*\).*/\1/'))

# Older versions of Go use the "-X foo bar" syntax.  Newer versions either warn
# or error on that syntax, and use "-X foo=bar" instead.
LINK_SEPERATOR=$(if $(filter 1.5, $(word 1, $(sort 1.5 $(GO_VERSION)))),=, )

.PHONY: pvc hmm-bc-cnv  report_single_partition len_mass_model_targeted len_mass_model lariat $(GOBINS) test clean cython hap.py

#
# Targets for development builds.
#
all: $(GOBINS) tada pvc hmm-bc-cnv report_single_partition len_mass_model len_mass_model_targeted lariat bamtofastq tenkit-all cython hap.py

marsoc: all version-files

version-files:
	git describe --tags --dirty --always > .version
	cd tenkit && git fetch --tags && git describe --tags --dirty --always > .version

tenkit-all:
	make -C tenkit all

lib/bin:
	mkdir lib/bin

lariat: lib/bin
	make -C $(LARIAT_DIR)/go lariat
	cp $(LARIAT_DIR)/go/bin/* lib/bin/

tada: lib/bin
	cd $(TADA_DIR); CARGO_HOME=$(TOP_DIR)/.cargo cargo build --release
	cp $(TADA_DIR)/target/release/tada lib/bin

$(GOBINS): lib/bin
	go install -ldflags "-X $@.__VERSION__$(LINK_SEPERATOR)'$(VERSION)'" $@
	cp lib/go/bin/* lib/bin

#
# Targets for cython builds
#
CYTHON_SRCS=$(shell find $(PWD)/mro/stages $(PWD)/lib/python -type f -name '*.pyx')
CYTHON_LIBS=$(patsubst %.pyx, %.so, $(CYTHON_SRCS))
CYTHON_BUILDPATH=$(BUILDPATH)/pkg/cython
CYTHON_FLAGS?=--line-directives $(EXTRA_CYTHON_FLAGS)

# Prevent make from automatically deleting intermediate files.
.PRECIOUS: $(CYTHON_BUILDPATH)/%.c $(CYTHON_BUILDPATH)/%.o

$(CYTHON_BUILDPATH)/%.c: $(PWD)/%.pyx
	mkdir -p $(@D) && cython $(CYTHON_FLAGS) -w $(<D) -o $(abspath $@) $(<F)

$(CYTHON_BUILDPATH)/%.o: $(CYTHON_BUILDPATH)/%.c
	$(CC) $(CFLAGS) -g -O3 -c -fPIC -fopenmp \
	    -I$(PYTHON_INC_DIR) \
	    -I$(PYTHON_LIB_DIR)/python2.7/site-packages/numpy/core/include \
	    -o $@ \
	    $<

$(PWD)/%.so: $(CYTHON_BUILDPATH)/%.o
	$(CC) $(LDFLAGS) -shared -fopenmp -L$(PYTHON_LIB_DIR) -lpython2.7 -fPIC $< -o $@

cython: $(CYTHON_LIBS)

clean-cython:
	rm -rf $(CYTHON_BUILDPATH)
	rm -f $(CYTHON_LIBS)

#
# Build STAN model binaries
#
STAN_DIR=lib/stan

$(STAN_DIR)/len_mass_model: $(STAN_DIR)/len_mass_model.stan
	$(MAKE) -C $(STAN_DIR)/cmdstan ../len_mass_model

$(STAN_DIR)/len_mass_model_targeted: $(STAN_DIR)/len_mass_model_targeted.stan
	$(MAKE) -C $(STAN_DIR)/cmdstan ../len_mass_model_targeted

lib/bin/len_mass_model: $(STAN_DIR)/len_mass_model lib/bin
	cp $< $@

lib/bin/len_mass_model_targeted: $(STAN_DIR)/len_mass_model_targeted lib/bin
	cp $< $@

len_mass_model: lib/bin/len_mass_model
len_mass_model_targeted: lib/bin/len_mass_model_targeted

clean-stan:
	rm -rf $(STAN_DIR)/len_mass_model $(STAN_DIR)/len_mass_model.cpp
	rm -rf $(STAN_DIR)/len_mass_model_targeted $(STAN_DIR)/len_mass_model_targeted.cpp
	find $(STAN_DIR) -name '*.o' | xargs rm -f
	find $(STAN_DIR) -name '*.a' | xargs rm -f

#
# Build PVC via rust cargo. 
# We set CARGO_HOME to .cargo in the longranger directory so 
# we can clean the carog cache to guarantee a fresh build.
#
TOP_DIR=$(shell pwd)
PVC_DIR=lib/pvc

pvc: lib/bin
	cd lib/pvc; CARGO_HOME=$(TOP_DIR)/.cargo cargo build --release
	cp lib/pvc/target/release/pvc lib/bin

clean-pvc:
	rm -Rf .cargo
	rm -Rf lib/bin/pvc
	rm -Rf lib/pvc/target

clean-tada:
	rm -Rf .cargo
	rm -Rf lib/bin/tada
	rm -Rf $(TADA_DIR)/target


# Build the rust version of report_single_partition 
# call it RRSP
TOP_DIR=$(shell pwd)

HBC_DIR=lib/rust/hmm-bc-cnv

hmm-bc-cnv: lib/bin
	cd ${HBC_DIR}; CARGO_HOME=$(TOP_DIR)/.cargo cargo build --release
	cp ${HBC_DIR}/target/release/hmm-bc-cnv lib/bin

clean-hmm-bc-cnv:
	rm -Rf .cargo
	rm -Rf lib/bin/hmm-bc-cnv
	rm -Rf ${HBC_DIR}/target

RRSP_DIR=lib/rust/report_single_partition

report_single_partition: lib/bin
	cd ${RRSP_DIR}; CARGO_HOME=$(TOP_DIR)/.cargo cargo build --release
	cp ${RRSP_DIR}/target/release/report_single_partition lib/bin

clean-report_single_partition:
	rm -Rf .cargo
	rm -Rf lib/bin/report_single_partition
	rm -Rf ${RRSP_DIR}/target

bamtofastq: lib/bin
	cd lib/bamtofastq; CARGO_HOME=$(TOP_DIR)/.cargo cargo build --release
	cp lib/bamtofastq/target/release/bamtofastq lib/bin

clean-bamtofastq:
	rm -Rf .cargo
	rm -Rf lib/bin/bamtofastq
	rm -Rf lib/bamtofastq/target

test: cython
	go test -v loupe/test
	MROPATH=$(PWD)/mro:$(PWD)/tenkit/mro \
		PATH=$(PATH):$(PWD)/bin \
		PYTHONPATH=$(PYTHONPATH):$(PWD)/lib/python \
		nosetests --with-xunit --with-coverage --cover-erase --cover-html --cover-xml --cover-package=stages,kitten mro/stages/* lib/python

clean-lariat:
	make -C $(LARIAT_DIR)/go clean

clean-post-tenkit-submodule:
	rm -rf $(LIBPY)/striped_smith_waterman
	rm -rf $(LIBPY)/tenkit
	rm -rf $(LIBPY)/tenx_numerics
	rm -rf mro/stages/bcl_processor
	rm -rf mro/stages/preflight/bcl_processor


