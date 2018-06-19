# Protein_DNA_Sequence_Tools
A repository for tools to analyze protein and/or DNA sequences

## Protein_Scrambler.py
A script that accepts a plain text file containing a single protein sequence and outputs a file with a scrambled version of the sequence.  
Necessary argument:  
The name of the file containing your protein sequence, in the same directory as this script or with the complete path name (-f, or --file).  
Optional arguments include;  
1) help (-h, or --help)
2) omitting the first residue from scrambling (-o, or --omit_start) if you don't wish to include the start Methionine in the scrambling  
3) specifying a subset of amino acids that you wish to scramble (-a, or --amino_acids)  
e.g. if you only want to scramble charged residues (D, E, H, K, and R), while maintaining other types of residues at their original positions.

#### Usage:
e.g. _python Protein_Scrambler.py [-h] -f sequence.txt [-o] [-a 'DEHKR']_
