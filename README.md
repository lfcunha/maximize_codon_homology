# maximize_codon_homology
Python script to maximize DNA sequence homology while retaining aminoacid sequence


it requires BioPython for reading inputs
Input:
    1) fasta file with multiple dna sequences
    2) .aln (clustalw) file with protein multiple seq alignment (with gaps)
    3) output prefix
    4) optional -v flag: prints result to console
This script maximizes sequence similarity of the 3rd nucleotide of each codon.
It optimizes not only same aminoacid codons, as well as codons coding different aa.
After optimziation, it validates that optimized dna sequence codes the same aa sequence as the 
original sequence, throwing an error otherwise. (It optimizes only the 3rd nt, not the first two) .

Usage: 
maximize_homology [-h] [--dna FASTA] [--protein_alignment ALIGNMENT] [--output OUTPUT_PREFIX]
