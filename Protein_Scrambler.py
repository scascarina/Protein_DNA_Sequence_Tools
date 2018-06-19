
import random
import sys
import argparse

def main():
    script = sys.argv[0]
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file',
                        help='Name of your file containing a single protein sequence')
                        
    parser.add_argument('-o','--omit_start', 
                        action='store_true',
                        help='Omit the first amino acid from scrambling (e.g. the start methionine)')
                        
    parser.add_argument('-a','--amino_acids',
                        default='ACDEFGHIKLMNPQRSTVWY',
                        help='Amino acids to include in the scrambling. When specified, only the indicated amino acids (e.g. "QNST") will be scrambled. Default includes all amino acids')
                        
    args = parser.parse_args()
    input_file = args.file
    
    h = open(input_file)
    sequence = ''
    for line in h:
        line = line.rstrip()
        sequence += line
    
    scrambled_seq = scrambler(sequence, args)
    
    output = open('scrambler_results.txt', 'w')
    output.write(scrambled_seq)
    output.close()
    
    
def scrambler(sequence, args):
    amino_acids = args.amino_acids
    
    #if omit_start is True, this ensures that the Methionine at the start of the protien isn't scrambled with the sequence
    if args.omit_start == True:
        print(sequence)
        sequence = sequence[1:]

    aa_to_randomize = ""            #initiates a list to store each occurrence of the user-specified amino acids to randomize
    for aa in sequence:             #For loop extracts specified amino acids from a protein sequence and creates a new string with them
        if aa in amino_acids:
            aa_to_randomize += aa
            
    # create a jumbled version of the extracted amino acids
    jumbled_extract = ""
    while aa_to_randomize:
        position = random.randrange(len(aa_to_randomize))
        jumbled_extract += aa_to_randomize[position]
        aa_to_randomize = aa_to_randomize[:position] + aa_to_randomize[(position + 1):]

    #Initializes seq2 variable. For loop builds a new string (assinged to seq2) that inserts the aa's from jumbled extract into the proper positions.
    seq2 = ''
    if args.omit_start == True:
        seq2 += 'M'
        
    for aa in sequence:
        if aa in amino_acids:
            seq2 += jumbled_extract[0]
            jumbled_extract = jumbled_extract[1:]
        else:
            seq2 += aa

    return seq2
    
if __name__ == '__main__':
    main()