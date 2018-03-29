#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#

.PHONY: all

all: target/release/tada

target/release/tada: src/*.rs src/*/*.rs
	RUSTFLAGS="${RUSTFLAGS} -A unused_must_use -A unused_variables -A unused_mut" cargo build --bin=tada --release


# BEGIN_INTERNAL - DO NOT MODIFY THIS LINE
# Everything above must be suitable for OSS release.
# Everything below will be stripped by oss-strip-makefile.

oss-strip-makefile:
	sed -n '/BEGIN_INTERNAL/q;p' Makefile > Makefile.oss

oss-ready: oss-strip-makefile 
	make -C external/rust-shardio oss-ready
	make -C external/rust-fastq-10x oss-ready
	make -C external/martian oss-ready
	mv Makefile.oss Makefile
	rm -rf .git
	rm -f  .gitmodules
	rm -f  .gitignore
	rm -f  .travis.yml
	rm -f  build.rs test_extern.py			# killing random stuff
	rm -f  README.md
	rm -rf play test_data target
	@echo repo oss ready

sake-strip-pd:
	rm -f  build.rs Cargo.lock Cargo.toml .gitignore Makefile README.md test_extern.py .travis.yml
	rm -rf external play src target test_data

sake-strip-cs: sake-strip-pd
	rm -f  mro/*_pd*.mro
	rm -rf mro/stages/*_pd*
	rm -rf mro/stages/**/*_pd*
