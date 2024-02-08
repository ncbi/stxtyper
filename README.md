# This is for development purposes only. Please do not use for research, public health, or god forbid diagnostic purposes because we are still developing and testing this software.

# StxTyper

StxTyper is used to determine stx type from nucleotide sequence. 

# Installation

## Compiling

```
git clone https://github.com/evolarjun/stxtyper.git
cd stxtyper
make
make test
```

# Usage

### Example

```
stxtyper -n nucleotide.fa
```

## Output

The output of StxTyper is a tab-delimited file with the following fields, all percent identity and coverage metrics are measured in proportion of amino-acids.

1. __target_contig__: The contig identifier from the input FASTA file
2. __stx_type__: The stx type called by the algorithm
3. __operon__: What status the operon was found to be. It can be
    - __STANDARD__ (soon to be "COMPLETE") for complete and fully typeable known stx types
    - __PARTIAL__ for partial operons that are internal to contigs and not terminating at contig boundaries
    - __PARTIAL_CONTIG_END__ for partial operons that could be split by contig boundaries due to sequencing or assembly artifacts
    - __INTERNAL_STOP__ for Stx operons where one of the subunits has a nonsense mutation
    - __FRAMESHIFT__ where StxTyper detected an indel in the coding sequence that would cause a frame shift in one or more of the subunits
    - __NOVEL__ a full-length stx operon that is not typeable using the current scheme
4. __identity__ The combined percent identity for both A and B subunits
5. __target_start__ The detected start of the alignments
6. __target_stop__ The detected end of the alignments
7. __target_strand__ What strand the target is on
8. __A_reference__ The closest reference protein for the A subunit, empty if none aligned
9. __A_identity__ The percent identity to the reference for the A subunit, empty if none aligned
10. __A_coverage__ The percentage of the reference for the A subunit that is covered by the alignment, empty if none aligned
11. _B_reference__ The closest reference protein for the B subunit, empty if none aligned
12. __B_identity__ The percent identity to the reference for the B subunit, empty if none aligned
13. __B_coverage__ The percentage of the reference for the B subunit that is covered by the alignment, empty if none aligned

