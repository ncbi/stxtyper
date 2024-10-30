# StxTyper methods description

## Only "complete" operons are typed

Only operons that have both stxA and stxB subunits aligning over 100% of their reference coding sequences are typed using the following algorithm.

## Typing for COMPLETE / NOVEL operons

- Compute the combined % identity as the sum of identities of the both proteins / sum of reference lengths of both proteins.
- A subtype is declared if
  - Both subunits are of the same stx type, considering stx2a, stx2c and stx2d as one generalized type “stx2acd”
  - Intergenic region between subunit A and subunit B is <= 36 bp
  - The combined % identity >= the cutoff which is:
    - 98.3 % for stx1 types
    - 98.5 % for stx2k and stx2l
   - 98 % for the other stx types
    - For the the types stx2a, stx2c, and stx2d, they are treated as one generalized type with a subtype declared if the operon has the following amino acids:

      ![image](https://github.com/user-attachments/assets/8252db28-6dbd-495c-b1e2-bdeaaad4f4a5)

      This table determines the exact stx type for the generalized “stx2acd” type, if an amino acid is present at those locations not on the chart then it may be called stx2, but is a novel subtype.
- If none of the above rules agree to define a subtype then it is a novel stx type. 

