from Bio import SeqIO

def make_reverse_complementary_sequence(seq):
    compBaseDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    compSeq = ''.join([compBaseDict[base] for base in seq])
    revCompSeq = compSeq[::-1]
    return revCompSeq

# NOTE: Porechop and filtlong occur here. It looks like they were not used in the first test. Circle back to later.
