
#dictionary with max side chain lengths in Angstroms
#defined as maximum distance from alpha carbon to furthest non-hydrogen atom in the sidechain
#taken from Table 1 in "Extent of N-terminal methionine excision from Esherichia coli proteins is governed by the side-chain length of the penultimate amino acid" Hirel et al (1989) PNAS
aa_lengths = {'A' : 1.51, 'C' : 2.83, 'D' : 3.74, 'E' : 4.97, 'F' : 5.10,
          'G' : 0.00, 'H' : 4.64, 'I' : 3.91, 'K' : 6.37, 'L': 3.90,
          'M' : 5.46, 'N' : 3.68, 'P' : 2.40, 'Q' : 4.93, 'R' : 7.40,
          'S' : 2.41, 'T' : 2.54, 'V' : 2.55, 'W' : 6.64, 'Y' : 6.43}

#dictionary with sidechain volumes (protein interior) in cubic angstroms
#defined as the volume of a given amino acid minus the volume of a glycine residue
#taken from column 2 in Table 3 from "Volume changes on protein folding" Harpaz et al (1994) Structure journal
#NOTE: these are aa sidechain volumes in protein interiors...all aa volumes are different in solution/on protein exterior
#NOTE: the value for cysteine assumes that it is not forming a disulfide bond...the volume if it has formed a disulfide bond is 39.7 cubic angstroms
aa_volumes = {'A' : 26.3, 'C' : 49.4, 'D' : 53.3, 'E' : 77.0, 'F' : 129.7,
          'G' : 0.0, 'H' : 95.5, 'I' : 101.1, 'K' : 106.2, 'L': 100.8,
          'M' : 103.9, 'N' : 63.7, 'P' : 59.3, 'Q' : 85.6, 'R' : 129.0,
          'S' : 30.4, 'T' : 56.2, 'V' : 75.3, 'W' : 167.9, 'Y' : 133.3}
              
def sidechain_LV():
    import sys
    from Bio import SeqIO
    script = sys.argv[0]
    file = sys.argv[1]
    
    output = open(file[:-6] + '_SideChainParameters_Results.csv', 'w')
    output.write('Amino Acid,Side-chain Length/Volume?,Values -->\n')
    
    #iterate through FASTA file records
    for seq_record in SeqIO.parse(file, 'fasta'):
        prot = str(seq_record.id)
        seq = str(seq_record.seq)

        #remove stop codon
        if seq.endswith('*'):
            seq = seq[:-1]

        #output residue-by-residue side-chain LENGTHS
        output.write(prot + ',Length')
        for aa in seq:
            if aa != '*':   #effectively skips internal stop codons...a few of the protein sequences in the yeast orf_trans.fasta file have internal stop codons.
                output.write(',' + str(aa_lengths[aa]))
        output.write('\n')
        
        #output residue-by-residue side-chain VOLUMES
        output.write(prot + ',Volume')
        for aa in seq:
            if aa != '*':   #effectively skips internal stop codons...a few of the protein sequences in the yeast orf_trans.fasta file have internal stop codons.
                output.write(',' + str(aa_volumes[aa]))
        output.write('\n')
        
    output.close()
    
    
if __name__ == '__main__':
    sidechain_LV()