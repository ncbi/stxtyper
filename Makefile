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
ifneq '$(CONDA_DB_DIR)' ''
	DBDIR := -D'CONDA_DB_DIR="$(CONDA_DB_DIR)"'
endif
ifneq '$(DEFAULT_DB_DIR)' ''
	DBDIR := -D'CONDA_DB_DIR="$(DEFAULT_DB_DIR)"'
endif

# for testing database updates using 
ifdef TEST_UPDATE
	TEST_UPDATE_DB := '-D TEST_UPDATE'
endif

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
COMPILE.cpp= $(CXX) $(CPPFLAGS) $(SVNREV) $(DBDIR) $(TEST_UPDATE_DB) -c 


.PHONY: all clean install release

BINARIES= stxtyper fasta_check db

all:	$(BINARIES)

release: clean
	svnversion . > version.txt
	make all

db: stx.prot
	makeblastdb  -in stx.prot  -dbtype prot > /dev/null

common.o:	common.hpp common.inc

stxtyper.o:  common.hpp common.inc 
stxtyperOBJS=stxtyper.o common.o tsv.o
stxtyper:	$(stxtyperOBJS)
	$(CXX) -o $@ $(stxtyperOBJS) -pthread $(DBDIR)

fasta_check.o:	common.hpp common.inc 
fasta_checkOBJS=fasta_check.o common.o 
fasta_check:	$(fasta_checkOBJS)
	$(CXX) -o $@ $(fasta_checkOBJS)


clean:
	rm -f *.o
	rm -f $(BINARIES)

install:
	@if [ ! -e $(DESTDIR)$(bindir) ]; \
	then \
		mkdir -p $(DESTDIR)$(bindir); \
	fi
	$(INSTALL) $(BINARIES) $(DESTDIR)$(bindir)

# stxtyper binaries for github binary release
GITHUB_FILE=stxtyper_binaries_v$(VERSION_STRING)
GITHUB_FILES = $(BINARIES)
github_binaries:
	@if [ ! -e version.txt ]; \
	then \
		echo >&2 "version.txt required to make a distribution file"; \
		false; \
	fi
#   first recompile stxtyper.o to pick up the new version info
	rm stxtyper.o stxtyper
	make
	mkdir $(GITHUB_FILE)
	echo $(VERSION_STRING) > $(GITHUB_FILE)/version.txt
	cp $(GITHUB_FILES) $(GITHUB_FILE)
	if [ -e $(GITHUB_FILE).tar.gz ]; then rm $(GITHUB_FILE).tar.gz; fi
	cd $(GITHUB_FILE); tar cvfz ../$(GITHUB_FILE).tar.gz *
#	tar cvfz $(GITHUB_FILE).tar.gz $(GITHUB_FILE)/*
	rm -r $(GITHUB_FILE)/*
	rmdir $(GITHUB_FILE)

#test : $(DISTFILES) Makefile *.cpp *.hpp *.inc test_dna.fa test_prot.fa test_prot.gff test_dna.fa test_dna.expected test_prot.expected test_both.expected
#	./test_stxtyper.sh
