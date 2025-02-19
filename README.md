# StxTyper

[![Latest StxTyper version](https://img.shields.io/github/v/release/ncbi/stxtyper)](https://github.com/ncbi/stxtyper/releases/latest) [![Latest Anaconda version](https://img.shields.io/conda/vn/bioconda/ncbi-stxtyper)](https://anaconda.org/bioconda/ncbi-stxtyper)

StxTyper is used to determine stx type from nucleotide sequence. Stx (Shiga-toxin) genes are found in some strains of _Escherichia coli_ and code for powerful toxins that can cause severe illness. StxTyper is software to classify these genes from assembled sequence using a standard algorithm.

## WARNING: This is currently beta software and changes and new releases may come quickly. Please report any issues or comments to pd-help@ncbi.nlm.nih.gov or <a href="https://github.com/ncbi/stxtyper/issues/new">open an issue</a> on GitHub.

# Installation

Note StxTyper is included with [AMRFinderPlus](https://github.com/ncbi/amr/wiki) as of version 4.0 and is run by AMRFinderPlus when the `--organism Escherichia` option is used. If you have installed AMRFinderPlus you don't need to separately install StxTyper.

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

## Docker

StxTyper is included with several Docker images including the [ncbi/amr](https://hub.docker.com/r/ncbi/amr/) docker image.

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

- `--amrfinder` Print the output in the fields that match AMRFinderPlus output. [See below for details](#--amrfinder-output).

- `--print_node` In the `--amrfinder` output format add the "Hierarchy node" as the last column. See the [field description in the AMRFinderPlus documentation](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#fields) for details.

- `--nucleotide_output <fasta_output_filename>` Print the nucleotide sequence of the identified operons to this file in FASTA format. Takes the sequence from _Target start_ to _Target stop_ and reverse-complements it if necessary to put it in the coding frame.

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

### `--amrfinder` output

This format of output matches the field names for AMRFinderPlus and is used when StxTyper is run as part of the AMRFinderPlus analysis pipeline. Note that AMRFinderPlus will include gene-based identification of Stx genes in separate rows. StxTyper will only output type calls for the operon as a whole.

1. __Protein id__ Always NA for StxTyper which runs using translated blast
   alignments against the nucleotide sequence.
2. __Contig id__ The FASTA identifier for the contig where the operon was
   found.
3. __Start__ The 1-based coordinate of the first base of the identified operon or partial operon.
4. __Stop__ The coordinate of the last base of the identified operon or partial operon.
5. __Strand__ '+' or '-' to indicate if the coding sequence of the genes is on the forward or reverse strand.
6. __Element symbol__ The operon symbol, corresponds to the "stx type" in the default output. For example 'stx2a_operon' or 'stx1_operon'
7. __Element name__ A description of the identified operon, including the subtype or whether it was a partial or frameshift operon.
8. __Scope__ Always 'plus' for StxTyper output corresponding to the [Scope](https://github.com/ncbi/amr/wiki/Interpreting-results#the-scope-column-and-the---plus-option) used for virulence genes by AMRFinderPlus.
9. __Type__ Always 'VIRULENCE' for StxTyper output.
10. __Subtype__ Always 'STX_TYPE' for StxTyper output.
11. __Class__ Either STX1 or STX2 corresponding to the type of the stx operon.
12. __Subclass__ The more detailed subtype of the operon if typeable.
13. __Method__ Corresponds to the __operon__ field in standard output. See above for details.
14. __Target length__ Calculated as __Stop__ - __Start__ or the length of the operon hit in nucleotide sequence.
15. __Referemce sequence length__ Always empty for StxTyper output.
16. __% Coverage of reference__ Always empty for StxTyper output.
17. __% Identity to reference__ The amino-acid percent identity to the reference genes, does not include the intergenic spacer.
18. __Alignment length__ The total amino-acid length of the subunit alignments.
19. __Closest reference accession__ The closest reference accessions, two values separated by a ', ' if both subunits aligned.
20. __Closest reference name__ The name of the closest reference operon.
21. __HMM accession__ Always NA for StxTyper output since HMMs are not used in operon typing.
22. __HMM description__ Always NA for StxTyper output
23. __Hierarchy node__ [Optional] When the `--print_node` option is used this is the nodes in the [Reference Gene Hierarchy](https://www.ncbi.nlm.nih.gov/pathogens/genehierarchy/) for each of the subunits separated by '::' if there are more than one (e.g., stxB2a::stxA2c). Note that for some Stx operon types the A and B subunits may have different types in the hierarchy becuase some subunits can appear in multiple stx types.

# Algorithm

## StxTyper operon detection
StxTyper (https://github.com/ncbi/stxtyper) uses translated BLAST (Camacho et al., 2009) alignments to a reference database of StxA and StxB protein sequences that is included with StxTyper and as part of the NCBI Pathogen Detection Reference Gene Catalog (https://www.ncbi.nlm.nih.gov/pathogens/refgene). Alignments are screened according to the following algorithm to determine whether they are full length, and only complete operons are typed down to the subtype (e.g., stx1c, stx2n) using the algorithm described below. Those that are not complete operons or cannot be fully typed by the scheme will be assigned stx types of stx1 or stx2, or if indeterminate, stx. 

StxTyper uses the following algorithm to identify putatively incomplete or non-functional operons. These categories are included in the “operon” field of StxTyper and also appear in the AMRFinderPlus “Method” field when run by AMRFinderPlus. Values of FRAME_SHIFT, INTERNAL_STOP, PARTIAL_CONTIG_END, EXTENDED, or PARTIAL indicate putatively incomplete or non-functional and therefore not subtyped operons. Operons that don’t meet any of these characteristics are designated either COMPLETE, COMPLETE_NOVEL, or AMBIGUOUS and, where possible given the complete operon typing scheme described below typed down to the subtype (e.g., stx2a). Stx operons that are putatively incomplete or non-functional are assigned to the type level (i.e., stx1 or stx2) if possible.

Algorithm for determining if operons are complete
Alignments or operons that are < 80% identical to individual reference subunits or < 80% identical to both A and B subunits combined are dropped from consideration. The putatively non-functional categories are determined in this order:

1.	FRAME_SHIFT is called when there are two consecutive BLAST hits to the same reference with distance < 10-bp in different reading frames. Because these blast matches to the whole reference protein are split into multiple alignments the numbers for percent identity and alignment length are calculated from the identified alignments and may not be what is expected.
2.	INTERNAL_STOP is called when there is a ‘*’ character (stop codon) any place in the query portion of the alignment other than at the last position of the reference sequence.
3.	PARTIAL_CONTIG_END is used to differentiate operons that may be split by assembly or sequencing issues at the end of a contig versus PARTIAL operons internal to contigs that are more likely to be incomplete in reality. PARTIAL_CONTIG_END operons are identified by alignments that terminate internal to the reference sequences < 3-bp from the end of a contig. A single-subunit operon is also identified as PARTIAL_CONTIG_END if the length of the unmatched portion of the contig in the same orientation as the missing subunit is < 36-bp (maximum allowed intergenic region) + 60-bp (minimum coding length to determine a subunit alignment has been missed) from the end of the contig. 
4.	EXTENDED is called when the whole reference protein aligns, but the stop codon is not present at the end of the alignment to the reference sequence.
5.	PARTIAL is called for the remaining cases where < 100% of the reference sequences of both subunits align.

## StxTyper typing algorithm

This typing scheme is only applicable to complete operons defined as 100% coverage over 100% of the reference protein sequence for both subunits that are < 36-bp apart. 

*	Compute the combined % identity as the sum of identities of the both proteins / sum of reference lengths of both proteins.
* A subtype is declared if
    * Both subunits are of the same stx type, considering stx2a, stx2c and stx2d as one generalized type “stx2acd”
    * Intergenic region between subunit A and subunit B is <= 36 bp
    * The combined % identity >= the cutoff which is:
        * 98.3 % for stx1 types
        * 98.5 % for stx2k and stx2l
        * 98 % for the other stx types
    * For the the types stx2a, stx2c, and stx2d, they are treated as one generalized type with a subtype declared if the operon has the following amino acids:

<table>
<tr><th>Stx type</th><th colspan=2>Subunit A</th><th>Subunit B</th></tr>
<tr><th>Position</th><th>313</th><th>319</th><th>35 (354 on holotoxin alignment)</th></tr>
<tr><th>stx2a</th><td>F/S</td><td>K/E</td><td>D</td></tr>
<tr><th>stx2c</th><td>F  </td><td>K/E</td><td>N</td></tr>
<tr><th>stx2d</th><td>S  </td><td>E  </td><td>N</td></tr>
</table>
    * This table determines the exact stx type for the generalized "stx2acd" type, if an amino acid is present at those locations not on the chart then it may be called stx2, but is a novel subtype.
* If none of the above rules agree to define a subtype then it is a novel stx type. 

