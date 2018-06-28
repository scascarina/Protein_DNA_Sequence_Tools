
"""
PAPA - a program for predicting the prion propensity of a protein

FoldIndex formula:
2.785(H) - |R| - 1.151
H = sum of the hydrophobicities across the window							
R = net charge (where D/E=-1; K/R=+1; all others = neutral (including H))							
"""

#window_size = 41
#half_window_size = int(window_size / 2)

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

propensities = {'A' :-0.396490246, 'C' : 0.415164505, 'D' : -1.276997939, 'E' : -0.605023827, 'F' : 0.838732498,
                'G' : -0.039220713, 'H': -0.278573356, 'I' : 0.813697862, 'K': -1.576748587, 'L' : -0.040005335,
                'M' : 0.673729095, 'N' : 0.080295334, 'P' : -1.197447496, 'Q' : 0.069168387, 'R' : -0.405858577,
                'S' : 0.133912418, 'T' : -0.11457038, 'V' : 0.813697862, 'W' : 0.666735081, 'Y' : 0.77865336}
																			
hydrophobicity = {'A' : 0.7, 'C' : 0.777777778, 'D' : 0.111111111, 'E' : 0.111111111, 'F' : 0.811111111,
                  'G' : 0.455555556, 'H' : 0.144444444, 'I' : 1, 'K' : 0.066666667, 'L': 0.922222222,
                  'M' : 0.711111111, 'N' : 0.111111111, 'P' : 0.322222222, 'Q' : 0.111111111, 'R' : 0,
                  'S' : 0.411111111, 'T' : 0.422222222, 'V' : 0.966666667, 'W' : 0.4, 'Y' : 0.355555556}

charge = {'A' : 0, 'C' : 0, 'D' : -1, 'E' : -1, 'F' : 0,
          'G' : 0, 'H' : 0, 'I' : 0, 'K' : 1, 'L': 0,
          'M' : 0, 'N' : 0, 'P' : 0, 'Q' : 0, 'R' : 1,
          'S' : 0, 'T' : 0, 'V' : 0, 'W' : 0, 'Y' : 0}


def window_scores(args, sequence, aa_dict, ignore_consecutive_prolines = False) :

    return [window_score(args, sequence, i, aa_dict, ignore_consecutive_prolines) for i in range(len(sequence))]

def window_score(args, sequence, position, aa_dict, ignore_consecutive_prolines = False) :
    #calculates average values for a given window, depending on which aa_dict is passed to this function
    #used to calculate average PAPA score, average hydrophobicity, and average charge
    start,end = get_window(args, sequence, position)
    score = 0.0
    for i in range(start, end) :
        if sequence[i] not in amino_acids :
            continue
        if (not ignore_consecutive_prolines) :
            score += aa_dict[sequence[i]]
        else :
            if (sequence[i] != 'P') :
                score += aa_dict[sequence[i]]
            elif ((i > 0 and sequence[i-1] == 'P') or (i > 1 and sequence[i-2] == 'P')) :
                pass
            else :
                score += aa_dict[sequence[i]]
    return score / (end - start)

def get_window(args, sequence, position) :
    start = max(position - args.half_window_size, 0)
    end = min(len(sequence), position + args.half_window_size + 1)
    return start,end
                                   

def fold_index(args, sequence) :
    charges = window_scores(args, sequence, charge)
    hydrophobicities = window_scores(args, sequence, hydrophobicity)
    fold_index_list = [2.785 * hydrophobicities[i]  - abs(charges[i]) - 1.151 for i in range(len(sequence))]
    return super_window_scores(args, sequence, fold_index_list)


def super_window_scores(args, sequence, window_scores, fold_index_scores = None) :
    scores = []
    for i in range(len(sequence) - args.window_size) :
        if (fold_index_scores is not None and fold_index_scores[i] > 0) :
            scores.append(None)
        else :
            score = 0.0
            weights = 0.0
            for j in range(i, i + args.window_size) :
                start,end = get_window(args, sequence, j)
                score += (end - start) * window_scores[j]
                weights += (end - start)
            scores.append(score / weights)
    return scores
    

def classify(args, sequence, ignore_fold_index = False) :

    fold_index_list = fold_index(args, sequence)
    #print 'number of positions with negative fold index', sum([1 for f in fold_index_list if f < 0]), 'out of', len(sequence)

    window_propensities = window_scores(args, sequence, propensities, True)
    if ignore_fold_index :
        scores = super_window_scores(args, sequence, window_propensities)
    else :
        scores = super_window_scores(args, sequence, window_propensities, fold_index_list)
        
    max_score = -1.0
    max_position = -1
    for score in scores:
        if score != None and score > max_score:
            max_score = score
            max_position = scores.index(max_score)

    return max_score, max_position, scores, fold_index_list
    

class MalformedInput :
    "Exception raised when the input file does not look like a fasta file."
    pass

class FastaRecord :
    "a fasta record."
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

def _fasta_itr(file_handle) :
    "Provide an iteration through the fasta records in file."

    h = file_handle.readline()[:-1]
    if h[0] != '>':
        raise MalformedInput()
    h = h[1:]

    seq = []
    for line in file_handle:
        line = line[:-1] # remove newline
        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))
            h = line[1:]
            seq = []
            continue
        seq.append(line)
    yield FastaRecord(h,''.join(seq))

class fasta_itr (object) :
    "An iterator through a sequence of fasta records."
    def __init__(self, src) :
        self.__itr = _fasta_itr(src)

    def __iter__(self) :
        return self

    def __next__(self) :
        return self.__itr.__next__()

#===================================================================================================================#

#Function Author: Sean Cascarina
#Added: 6/27/2018

def merge_overlapping(hit_positions):
    """A function that accepts a list of positions that are above a certain PAPA threshold and
        effectively merges positions that would result in overlapping PAPA-positive windows.
        
        The return is a new list of 2-tuples, where the first value is the start of the PAPA-positive window, 
        and the second value is the end of the PAPA-positive window
        e.g. if the original position list = [1, 2, 3, 4, 5, 8, 75, 115, 116, 117, 234, 236, 237],
        then this function would return the list [(0,177) , (213,297)]
        
        The position_extension_combinations list is an intermediate list that is used to calculate
        the start and end positions of the full PAPA-positive windows.
        This is a list of 2-tuples, where the first value is the position, and the second value
        is the number of extension residues to add to the position.
        The position_extension_combinations list
        e.g. for a hypothetical 2-tuple (136,5), the full PAPA window would be calculated as amino acids 116-222.
        The calculation would be (136-20) + 81 + 5
                                     ^      ^    ^
                                     |      |    |
                           start position   |   extension residues
                                      papa window size
                                      
        """
    
    j = 0
    position_extension_combinations = []
    
    #Runs until hit_positions is empty. Every time a new merged window is generated, that window is removed from hit_positions.
    #This eventually results in an empty hit_positions list.
    while hit_positions:
        start = hit_positions[0]
        
        #adds one to the index j while there are at least two items in the list and the two values at j and j+1 are within 81 amino acids of each other, indicating overlapping PAPA windows
        #once these conditions aren't met, the loop exits and appends a tuple to position_extension_combinations with the start position as the first value and the number of extension residues as the second value (extension residues are the number of amino acids to add to the start position to get to the position of the last score with an overlapping window)
        while j < (len(hit_positions) - 1) and ( (hit_positions[j+1] - 21) <= (hit_positions[j] + 60) ):
            j += 1
        
        # Gets the ending position for the current window and outputs the 2-tuple.
        end = hit_positions[j]
        extension_residues = end - start
        position_extension_combinations.append( (start , extension_residues) )

        #modifies hit_positions to delete all positions that were just merged, then re-sets j=0 to start iterating at the first position in the new list
        hit_positions = hit_positions[j+1:]
        j = 0
        
    return position_extension_combinations

#================================================================================================================#    

#Function Author: Sean Cascarina
#Added: 6/27/2018
    
def get_high_scoring_indices(sequence, highest_score_position, papa_scores, threshold=0.05):

    hit_indices = []
    
    #get window indices for all prlds above threshold
    filtered_nones = [float(x) for x in papa_scores if x != None ]
    if len(filtered_nones) > 0 and float(max( filtered_nones )) >= threshold:
        index = 0
        for score in papa_scores:
            if score != None and float(score) >= threshold:
                hit_indices.append(index)
            index += 1
        position_extension_combinations = merge_overlapping(hit_indices)
        
        #use position/extensions combinations to calculate the actual indices of the start and end of PAPA positive windows
        hit_indices = [ ((x-21) , (x+60+y)) if x>=21 else (0 , (x+60+y)) for (x,y) in position_extension_combinations ] #corrects for positions that are near the N-terminus to prevent negative x indices
        hit_indices = [ (x,y) if y <= len(sequence[:-1]) else ( x, len(sequence[:-1]) ) for (x,y) in hit_indices ] #corrects for positions that are near the C-terminus to prevent y indices that extend beyond the length of the protein. NOTE: the sequences in orf_trans.dat still contain the '*' stop codon.

    #otherwise, output the window indices for the highest scoring 81aa region only
    else:
        position_extension_combinations = [ (int(highest_score_position), 0) ]
        
        #use position/extensions combinations to calculate the actual indices of the start and end of PAPA positive windows
        hit_indices = [ ((x-21) , (x+60+y)) if x>=21 else (0 , (x+60+y)) for (x,y) in position_extension_combinations ] #corrects for positions that are near the N-terminus to prevent negative x indices
        hit_indices = [ (x,y) if y <= len(sequence[:-1]) else ( x, len(sequence[:-1]) ) for (x,y) in hit_indices ] #corrects for positions that are near the C-terminus to prevent y indices that extend beyond the length of the protein. NOTE: the sequences in orf_trans.dat still contain the '*' stop codon.

    return hit_indices
    
#=========================================================================================================================#

#Function Author: Sean Cascarina
#Added: 6/27/2018

def get_merged_window_seqs(sequence, hit_indices):
    """A function that accepts the name of the gene of interest and a list of the indices of all PAPA-positive regions
        and returns the amino acid sequences of those regions"""
        
    seqs = []
    for pos in hit_indices:
        start, end = pos
        seqs.append(sequence[ start : end] )

    return seqs
    
#=====================================================================================================#

#Modified by: Sean Cascarina
#Modification Date: 6/27/2018

def run(args) :

    if args.outfile is None :
        outfile = sys.stdout
    else :
        outfile = open(args.outfile, 'w')
        
    if args.window_size != 41:
        print('\n\nThis script only supports the default 41 amino acid window size.\n\n')
        exit()

    outfile.write('FASTA Header\tHighest PAPA Score\tPosition of Highest PAPA Score\tSequences of Highest-scoring Window(s)\tBoundaries of Highest-scoring Window(s)\n')
    for fasta_record in fasta_itr(open(args.fasta_input)) :
        sequence_id = str(fasta_record.header)
        sequence = str(fasta_record.sequence.upper())
        if len(sequence) <= args.window_size :
            outfile.write(sequence_id + '\t' + 'protein length below window size\n')
            continue
        
        high_score,pos,scores,fold_indexes = classify(args, sequence, args.ignore_fold_index)
        
        if high_score == -1.0:
            outfile.write(sequence_id + '\t' + str(high_score) + '\n')
            continue   
            
        #get the window indices for all PAPA-positive windows and the corresponding sequences
        window_indices = get_high_scoring_indices(sequence, pos, scores)
        hit_sequences = get_merged_window_seqs(sequence, window_indices)
        
        #shorten papa_score values by rounding floats to 5 places past the decimal, and by shortening 'None' to 'N'
        #shortening the values is necessary because the papa_scores list for very long proteins exceeds the maximum allowance in excel files
        #NOTE: this works for all yeast proteins, but may be an issue in other organisms with longer proteins
        papa_scores = [ round(float(x), 5) if x != None else 'N' for x in scores ]
        
        #START OUTPUT===========================
        outfile.write(sequence_id + '\t' + str(high_score) + '\t' + str(pos) + '\t')
        
        #Joins the PAPA-positive (or maximum score) sequences with a _;_ delimiter (underscore/semicolon/underscore) for visual clarity and easy parsing
        #and output sequences
        if len(hit_sequences) > 1:
            seqs_token = ['_;_'.join(seq for seq in hit_sequences)]
            outfile.write(str(seqs_token)[2:-2] + '\t')
        else:
            outfile.write(str(hit_sequences)[2:-2] + '\t')
            
        #Joins the PAPA-positive (or maximum score) window indices with a _;_ delimiter (underscore/semicolon/underscore) for visual clarity and easy parsing
        #and output window indices
        if len(window_indices) > 1:
            windows_token = ['_;_'.join(str( (window[0]+1, window[1]+1) ) for window in window_indices)]
            outfile.write(str(windows_token)[2:-2] + '\n')
        else:
            windows_token = [(window_indices[0][0]+1 , window_indices[0][1]+1)]
            outfile.write(str(windows_token)[1:-1] + '\n')

    
def test() :
    scores = {'Ure2' : 0.1031, 'Sup35': 0.0997, 'Rnq1' : 0.1413, 'New1' : 0.1398, 'Nsp1' : 0.0239, 'Puf2' : 0.0768, 'Pub1' : 0.1551}
    sequences = {}
    for fasta_record in fasta_itr(open('sequences.fasta')) :
        sequences[fasta_record.header] = fasta_record.sequence
    for id in sequences :
        score,pos,scores,fold_index = classify(sequences[id])
        print(id, len(sequences[id]), pos, score, scores[id], score / scores[id])
    print('scores ignoring foldIndex')
    for id in sequences :
        score,pos,scores,fold_index = classify(sequences[id], True)
        print(id, len(sequences[id]), pos, score, scores[id], score / scores[id])

#==================================================================================================================#

#Modified by: Sean Cascarina
#Modification Date: 6/28/2018

def parse_arguments(arguments) :
    import argparse
    parser = argparse.ArgumentParser(description='Predict whether a given protein is prion forming', prog = 'papa')
    parser.add_argument('fasta_input', help = 'the input file (in fasta format)')
    parser.add_argument('-o', '--outfile', type = str,
                        help = """the output file. Columns in the tab-delimited output file are:
                        1) Sequence id, 
                        2) Maximum score, 
                        3) Position of maximum score
                        4) Sequence(s) of Highest-scoring Window(s),
                        5) Boundaries of Highest-scoring Window(s)
                        
                        If the maximum score is below the default PAPA threshold 0.05, the sequence and boundaries
                        of the highest-scoring window will be the 81aa window corresponding to the
                        position of maximum score
                        
                        NOTE: this script does not support a "verbose" mode. Rather, it uses the data 
                        derived from verbose mode of the original PAPA script to define the 
                        boundaries of all high-scoring protein regions.
                        
                        This script also only supports the default 41 amino acid window size.
                        """)
    #parser.add_argument('-v', '--verbose', action='store_true', help = 'in verbose mode the output includes the scores for every position as well as the per-position fold index scores')
    parser.add_argument('--ignore_fold_index', action='store_true',
                        help = """Whether to ignore the fold index values for the protein (default: False)
                        Set it to True to ignore fold index.
                        """)
    parser.add_argument('--window_size', type = int, default = 41,
                        help = """The window size used by PAPA (default = 41).
                        Proteins that are shorter than the window size are not classified.
                        """)
    args = parser.parse_args(arguments)
    args.half_window_size = int(args.window_size / 2)    
    return args

if __name__ == '__main__' :

    import sys
    args = parse_arguments(sys.argv[1:])
    run(args)
    #test()
