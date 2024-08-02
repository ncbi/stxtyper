##############################################################################
# PUBLIC DOMAIN NOTICE This software/database is "United States Government
# Work" under the terms of the United States Copyright Act. It was written as
# part of the authors' official duties for the United States Government and
# thus cannot be copyrighted. This software/database is freely available to the
# public for use without a copyright notice. Restrictions cannot be placed on
# its present or future use.
#
# Although all reasonable efforts have been taken to ensure the accuracy and
# reliability of the software and data, the National Center for Biotechnology
# Information (NCBI) and the U.S. Government do not and cannot warrant the
# performance or results that may be obtained by using this software or data.
# NCBI, NLM, and the U.S. Government disclaim all warranties as to performance,
# merchantability or fitness for any particular purpose.
#
# In any work or product derived from this material, proper attribution of the
# authors as the source of the software or data should be made, using:
# https://ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/ as the
# citation.
###############################################################################

# the SVNREV is set automatically here for convenience,
# but when actually building we should override it like:
# make all SVNREV=-D\'SVN_REV=\"$VERSION\"\' or use
# a version.txt file
ifeq ($(wildcard version.txt),)
	VERSION_STRING := $(shell git describe --tags)
else
	VERSION_STRING := $(shell cat version.txt)
endif
SVNREV := -D'SVN_REV="$(VERSION_STRING)"'

INSTALL=install

# make it possible to hard define a database directory
# Define default paths
# This is a little convoluted because I broke things and don't want
# to change two different ways of defining the paths. This could
# be simplified in a later release
PREFIX ?= /usr/local
ifneq '$(INSTALL_DIR)' ''
	bindir=$(INSTALL_DIR)
endif
bindir ?= $(PREFIX)/bin

# detect system architecture and set appropriate flags
# this is probably not the best way (i.e. M1 Mac would be arm64)
# but it works for Nvidia Jetson boards (aarch64) 
ARCH := $(shell uname -m)
OS := $(shell uname -s)
# "hack": if amd64 we can set to aarch64
# as AArch64 and ARM64 refer to the same thing
# this should build for Mac M1 and other arm64 chips
ifeq ($(ARCH),arm64)
  ARCH := aarch64
endif
# report detected OS and arch in stdout
$(info Detected architecture: $(OS) $(ARCH))
# set CFLAGS based on arch
ifeq ($(ARCH),aarch64)
  # set arm CFLAGS
  CPPFLAGS = -std=gnu++17 -pthread --signed-char -falign-jumps -fno-math-errno -O3 
else
  # set x86_x64 CFLAGS
  CPPFLAGS = -std=gnu++17 -pthread -malign-double -fno-math-errno -O3
endif
# was: -std=gnu++14 

CXX=g++
COMPILE.cpp= $(CXX) $(CPPFLAGS) $(SVNREV) $(DBDIR) -c 


.PHONY: all clean install release test

BINARIES= stxtyper fasta_check fasta_extract
DATABASE= stx.prot
TEST= test/amrfinder_integration.expected   test/cases.expected \
test/amrfinder_integration.fa         test/cases.fa \
test/amrfinder_integration2.expected  test/synthetics.expected \
test/amrfinder_integration2.fa        test/synthetics.fa \
test/basic.expected                   test/virulence_ecoli.expected \
test/basic.fa                         test/virulence_ecoli.fa \
test/basic.nuc_out.expected


all:	$(BINARIES)

#db: stx.prot
#	makeblastdb  -in stx.prot  -dbtype prot > /dev/null

common.o:	common.hpp common.inc

stxtyper.o:  common.hpp common.inc 
stxtyperOBJS=stxtyper.o common.o tsv.o
stxtyper:	$(stxtyperOBJS)
	$(CXX) -o $@ $(stxtyperOBJS) -pthread $(DBDIR)

fasta_check.o:	common.hpp common.inc 
fasta_checkOBJS=fasta_check.o common.o 
fasta_check:	$(fasta_checkOBJS)
	$(CXX) -o $@ $(fasta_checkOBJS)

fasta_extract.o:	common.hpp common.inc 
fasta_extractOBJS=fasta_extract.o common.o 
fasta_extract:	$(fasta_extractOBJS)
	$(CXX) -o $@ $(fasta_extractOBJS)

clean:
	rm -f *.o
	rm -f $(BINARIES)

install:
	@if [ ! -e $(DESTDIR)$(bindir) ]; \
	then \
		mkdir -p $(DESTDIR)$(bindir); \
	fi
	$(INSTALL) $(BINARIES) $(DATABASE) $(DESTDIR)$(bindir) 

# stxtyper binaries for github binary release
GITHUB_FILE=stxtyper_v$(VERSION_STRING)_$(ARCH)_$(OS)
GITHUB_FILES = $(BINARIES) $(DATABASE) version.txt test_stxtyper.sh
binary_dist:
	@if [ ! -e version.txt ]; \
	then \
		echo >&2 "version.txt required to make a distribution file"; \
		false; \
	fi
	mkdir $(GITHUB_FILE)
	cp $(GITHUB_FILES) $(GITHUB_FILE)
	mkdir $(GITHUB_FILE)/test
	cp $(TEST) $(GITHUB_FILE)/test
	if [ -e $(GITHUB_FILE).tar.gz ]; then rm $(GITHUB_FILE).tar.gz; fi
	tar cvfz $(GITHUB_FILE).tar.gz $(GITHUB_FILE)/*
	rm -r $(GITHUB_FILE)/*
	rmdir $(GITHUB_FILE)

test : $(DISTFILES) Makefile *.cpp *.hpp *.inc test/* $(BINARIES)
	./test_stxtyper.sh
