# StxTyper

StxTyper is used to determine stx type from nucleotide sequence. Stx (Shiga-toxin) genes are found in some strains of _Escherichia coli_ and code for powerful toxins that can cause severe illness. StxTyper is software to classify these genes from assembled sequence using a standard algorithm.

## WARNING: This is currently beta software and changes and new releases may come quickly. Please report any issues or comments to pd-help@ncbi.nlm.nih.gov or open an issue on GitHub.

# Installation

## Installing with Bioconda

You'll need Mamba ([Installation instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)) first.

    micromamba create -n stxtyper ncbi-stxtyper \
      --channel conda-forge \
      --channel bioconda \
      --channel defaults \
      --strict-channel-priority

## Installing from binary

### Prerequisites 

#### NCBI BLAST+

See below under "Compiling" for instructions to install NCBI BLAST+.

### Download and install Binary

Download the latest binary tarball from https://github.com/ncbi/stxtyper/releases. Untar it and run the tests e.g.:

    tar xvfz stxtyper_v*.tar.gz
    cd stxtyper_v*/
    ./test_stxtyper.sh

Note that we are currently only publishing binary tarballs for x86 linux.

## Compiling

### Prerequisites

#### NCBI BLAST+

StxTyper needs NCBI BLAST binaries in your path (specifically tblastn). If you
don't already have BLAST installed see
https://www.ncbi.nlm.nih.gov/books/NBK569861/ for the official instructions to
install BLAST binaries. It's also available in many package repositories, for
example on Ubuntu:

    sudo apt-get install ncbi-blast+

#### C compiler and make

These are necessary if compiling from source. If using the binary distribution,
or Bioconda you won't need to worry about these. They generally come standard
for unix systems, if not the user will need to intall make and GCC. MacOS users
will need to go to the [App store and install
Xcode](https://apps.apple.com/in/app/xcode/id497799835?mt=12). 

### Compiling

StxTyper should compile cleanly for Mac and Linux x86 and ARM, though our official policy is we only support x86 Linux.

    git clone https://github.com/evolarjun/stxtyper.git
    cd stxtyper
    make
    make test


# Usage

    stxtyper -n <assembled_nucleotide.fa> [<options>]

### Example

    stxtyper -n nucleotide.fa

## Parameters

- `-nucleotide <nucleotide_fasta>` or `-n <nucleotide_fasta>` Assembled nucleotide sequence to search in FASTA format.

- `--name <assembly_identifier>` Add an identifier as the first column in each row of the report. This is useful when combining results for many assemblies.

- `--output <output_file>` or `-o <output_file>` Write the output to \<output\_file\> instead of STDOUT

- `--nucleotide_output <output_fasta>` Output the nucleotide sequence for any identified stx operons (includes partial and full length operons)

- `--blast_bin <path>` Directory to search for tblastn binary. Overrides environment variable `$BLAST_BIN` and the default PATH.

- `-q` or `--quiet` Suppress the status messages normally written to STDERR.

- `--log <log_file>` Error log file, appended and opened when you first run the application. This is used for debugging.

- `--debug` Run in debug mode. Additional messages are printed and files in $TMPDIR are not removed after running.

### For AMRFinderPlus

These options are not expected to be used outside of the [AMRFinderPlus](https://github.com/ncbi/amr/wiki) pipeline.

- `--amrfinder` Output in [AMRFinderPlus](https://github.com/ncbi/amr/wiki) format

- `--print_node` Add column for [Hierarchy node](https://www.ncbi.nlm.nih.gov/pathogens/docs/gene_hierarchy/) optionally reported by AMRFinderPlus.

## Output

The output of StxTyper is a tab-delimited file with the following fields, all percent identity and coverage metrics are measured in proportion of amino-acids.

1. __target_contig__: The contig identifier from the input FASTA file
2. __stx_type__: The stx type called by the algorithm, for "operon = COMPLETE"
   it will be stx plus two characters (e.g., stx1a), for other values of operon
   stx_type will be stx1, stx2, or just stx if it can't resolve at all.
3. __operon__: What status the operon was found to be. It can be
    - __COMPLETE__ for complete and fully typeable known stx types
    - __PARTIAL__ for partial operons that are internal to contigs and not
      terminating at contig boundaries
    - __PARTIAL_CONTIG_END__ for partial operons that could be split by contig
      boundaries due to sequencing or assembly artifacts
    - __EXTENDED__ The coding sequence extends beyond the reference stop codon
      for one or both of the reference proteins
    - __INTERNAL_STOP__ for Stx operons where one of the subunits has a
      nonsense mutation
    - __FRAMESHIFT__ where StxTyper detected an indel in the coding sequence
      that would cause a frame shift in one or more of the subunits
    - __AMBIGUOUS__ StxTyper found an ambiguous base in the query sequence
      (e.g., N), this could be the result sequencing or assembly error so the
      user might want to take a closer look at the sequence.
    - __COMPLETE_NOVEL__ a full-length stx operon that is not typeable using
      the current scheme
4. __identity__ The combined percent identity for both A and B subunits
5. __target_start__ The detected start of the alignments
6. __target_stop__ The detected end of the alignments
7. __target_strand__ What strand the target is on
8. __A_reference__ The closest reference protein for the A subunit, empty if
none aligned
9. __A_identity__ The percent identity to the reference for the A subunit,
empty if none aligned
10. __A_reference_subtype__ The subtype assigned to the reference sequence for
the A subunit. Note this may be different from the subtype for the operon as a
whole.
11. __A_coverage__ The percentage of the reference for the A subunit that is
covered by the alignment, empty if none aligned
12. __B_reference__ The closest reference protein for the B subunit, empty if
none aligned
13. __B_reference_subtype__ The subtype assigned to the reference sequence for
the B subunit. Note this may be different from the subtype for the operon as a
whole.
14. __B_identity__ The percent identity to the reference for the B subunit,
empty if none aligned
15. __B_coverage__ The percentage of the reference for the B subunit that is
covered by the alignment, empty if none aligned

