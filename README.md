# This is for development purposes only. Please do not use for research, public health, or god forbid diagnostic purposes because we are still developing and testing this software.

# StxTyper

StxTyper is used to determine stx type from nucleotide sequence. 

# Installation

## Prerequisites

### C compiler and make

These generally come standard for unix systems, if not the user will need make
and GCC. MacOS users will need to go to the [App store and install
Xcode](https://apps.apple.com/in/app/xcode/id497799835?mt=12). 

### NCBI BLAST

StxTyper needs NCBI BLAST binaries in your path (specifically tblastn). If you don't
already have BLAST installed see https://www.ncbi.nlm.nih.gov/books/NBK569861/
for instructions to install BLAST binaries.

Blast can also be installed using bioconda by first installing bioconda, then making sure to activate the environment in which it's installed. 
- See the [official install miniconda instructions for Linux](https://docs.anaconda.com/free/miniconda/#quick-command-line-install)
- See the [official install miniconda instructions for MacOS](https://docs.conda.io/en/latest/miniconda.html)

Then run:

    source ~/miniconda3/bin/activate
    conda create -y -c conda-forge -c bioconda -n blast blast
    conda activate blast

If you install BLAST via conda in this way you will need to run `conda activate blast` before you can run StxTyper.

## Compiling

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

- `--blast_bin <path>` Directory to search for tblastn binary. Overrides environment variable `$BLAST_BIN` and the default PATH.

- `-q` or `--quiet` Suppress the status messages normally written to STDERR.

- `--log <log_file>` Error log file, appended and opened when you first run the application. This is used for debugging

## Output

The output of StxTyper is a tab-delimited file with the following fields, all percent identity and coverage metrics are measured in proportion of amino-acids.

1. __target_contig__: The contig identifier from the input FASTA file
2. __stx_type__: The stx type called by the algorithm, for "operon = COMPLETE"
   it will be stx plus two characters (e.g., stx1a), for other values of operon stx_type will be
   stx1, stx2, or just stx if it can't resolve at all.
3. __operon__: What status the operon was found to be. It can be
    - __COMPLETE__ for complete and fully typeable
      known stx types
    - __PARTIAL__ for partial operons that are internal to contigs and not
      terminating at contig boundaries
    - __PARTIAL_CONTIG_END__ for partial operons that could be split by contig
      boundaries due to sequencing or assembly artifacts
    - __EXTENDED__ Where the coding sequence would extend beyond the stop codon
      for the reference protein
    - __INTERNAL_STOP__ for Stx operons where one of the subunits has a
      nonsense mutation
    - __FRAMESHIFT__ where StxTyper detected an indel in the coding sequence
      that would cause a frame shift in one or more of the subunits
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
10. __A_coverage__ The percentage of the reference for the A subunit that is
covered by the alignment, empty if none aligned
11. __B_reference__ The closest reference protein for the B subunit, empty if
none aligned
12. __B_identity__ The percent identity to the reference for the B subunit,
empty if none aligned
13. __B_coverage__ The percentage of the reference for the B subunit that is
covered by the alignment, empty if none aligned

