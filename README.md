# Sequence_Tools
A repository for tools to analyze protein and/or DNA sequences
<br>
<br>

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
e.g. _python Protein_Scrambler.py [-h] -f "sequence.txt" [-o] [-a 'DEHKR']_
<br>
<br>
<br>

## Sidechain_Length_Volume.py
A script that accepts a FASTA file containing any number of protein sequences (in FASTA format) and outputs two lines per protein. The first line contains residue-by-residue side-chain length data, whereas the second line contains residue-by-residue side-chain volume data.

Necessary argument:  
The name of the file containing your protein sequence(s), in the same directory as this script or with the complete path name (-f, or --file).  

#### Usage:
e.g. _python Sidechain_Length_Volume.py [-h] -f "sequences.fasta"_

##### References:
Side chain lengths were taken from Table 1 in: "Extent of N-terminal methionine excision from Esherichia coli proteins is governed by the side-chain length of the penultimate amino acid" Hirel et al (1989) _PNAS_  

Side chain volumes were taken from Table 3 (column 3; SCp) in: "Volume changes on protein folding" Harpaz et al (1994) _Structure_
<br>
<br>
<br>

## papa_Python3.py
A script that accepts a FASTA file containing any number of protein sequences (in FASTA format) and scores the prion propensity of each protein. This script was originally written by [Asa Ben-Hur](https://www.cs.colostate.edu/~asa/), but modified by me for compatibility with Python3. I also modified the output to be a tab-delimited file for simpler parsing. The output file can be opened with Excel and should be automatically parsed (although this may depend on Excel version/operating system) to separate columns for the proteinID, PAPA score, and location of the highest-scoring window.

#### Usage:
Refer to the [PAPA online server documentation](http://combi.cs.colostate.edu/supplements/papa/readme.html) for usage.
You can also use the [PAPA online server](http://combi.cs.colostate.edu/supplements/papa/) to analyze individual protein sequences.

##### References:
If you use PAPA, please cite: **Compositional determinants of prion formation in yeast**, Toombs _et al._ 2010 MCB [pubmed link](https://www.ncbi.nlm.nih.gov/pubmed/19884345)
