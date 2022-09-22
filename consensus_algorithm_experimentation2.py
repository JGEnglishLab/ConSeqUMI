from Bio import SeqIO, pairwise2
from statistics import median
import numpy as np
import itertools
from collections import Counter
from difflib import SequenceMatcher
from functools import reduce # Python 3

def most_common(lst):
    data = Counter(lst)
    return max(lst, key=data.get)

def initialize_consensus_output_string(strs, excerpt_length, start=True):
    if start: tempList = [x[:excerpt_length] for x in strs]
    else: tempList = [x[-excerpt_length:] for x in strs]
    final = max(set(tempList), key = tempList.count)
    return final

def overlap(a, b):
    temp = [i for i in range(1, len(b)+1) if a.endswith(b[:i])]
    if temp: return max(temp)
    return len(b)-1

class ConsensusSequenceGenerator():

    def __init__(self):
        self.binSize = 14
        self.consLength = 10

    def generate_consensus_sequence(self, file):

        seqStrs = [str(record.seq) for record in SeqIO.parse(file, "fastq")]
        self.numBins = round(median(map(len, seqStrs))/self.binSize)
        consWindows = []
        for seq in seqStrs:
            windows = [seq[i:i + self.consLength] for i in range(len(seq) - self.consLength + 1)]
            consWindows.append(np.array(windows))

        windows = [self.find_consensus_windows(consWindows, offset = i) for i in range(0, 15, 3)]
        allSeqs = []
        for i in range(len(windows[0])-1):
            seqs = [x[i] for x in windows]
            seqs.append(windows[0][i+1])
            tempSeq = reduce(lambda a, b: a + b[overlap(a,b):], seqs)
            #if not tempSeq.endswith(windows[0][i+1]): return 'ERROR'
            allSeqs.append(tempSeq)
        consensusSeq = reduce(lambda a, b: a + b[overlap(a,b):], allSeqs)
        return consensusSeq



    def return_possible_consensus_indices(self, l, offset):
        """Yield n number of sequential chunks from l."""
        d, r = divmod(len(l), self.numBins)
        numBinsRange = range(self.numBins)
        if offset: numBinsRange = range(self.numBins - 1)
        for i in numBinsRange:
            si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
            diff = (d+1 if i < r else d)
            yield np.array(range(si, si + diff - self.consLength + 1)) + offset


    def find_consensus_windows(self, seqWindows, offset = 0):
        binnedSeqs = []
        for window in seqWindows:
            if len(window) / self.numBins < self.consLength: continue
            tempSets = [set(window[x]) for x in self.return_possible_consensus_indices(window, offset)]
            binnedSeqs.append(tempSets)

        consSeqs = []

        numBinsRange = range(self.numBins)
        if offset: numBinsRange = range(self.numBins - 1)
        for binIdx in numBinsRange:
            allPossibleCons = [list(x[binIdx]) for x in binnedSeqs]
            allPossibleCons = list(itertools.chain.from_iterable(allPossibleCons))
            consSeqs.append(most_common(allPossibleCons))
        #print('--'.join(consSeqs))
        return consSeqs
'''
binFile = 'scratchpaper_code/seq_bin0.fq'
csg = ConsensusSequenceGenerator()
consSeq = csg.generate_consensus_sequence(binFile)
cons1 = 'TGCCACCCGGCTTCAACGAGTACGACTTCGTGCCCGAGAGCTTCGACCGGGACAAAACCATCGCCCTGATCATGAACAGTAGTGGCAGTACCGGATTGCCCAAGGGCGTAGCCCTACCGCACCGCACCGCTTGTGTCCGATTCAGTCATGCCCGCGACCCCATCTTCGGCAACCAGATCATCCCCGACACCGCTATCCTCAGCGTGGTGCCATTTCACCACGGCTTCGGCATGTTCACCACGCTGGGCTACTTGATCTGCGGCTTTCGGGTCGTGCTCATGTACCGCTTCGAGGAGGAGCTATTCTTGCGCAGCTTGCAAGACTATAAGATTCAATCTGCCCTGCTGGTGCCCACACTATTTAGCTTCTTCGCTAAGAGCACTCTCATCGACAAGTACGACCTAAGCAACTTGCACGAGATCGCCAGCGGCGGGGCGCCGCTCAGCAAGGAGGTAGGTGAGGCCGTGGCCAAACGCTTCCACCTACCAGGCATCCGCCAGGGCTACGGCCTGACAGAAACAACCAGCGCCATTCTGATCACCCCCGAAGGGGACGACAAGCCTGGCGCAGTAGGCAAGGTGGTGCCCTTCTTCGAGGCTAAGGTGGTGGACTTGGACACCGGCAAGACACTGGGTGTGAACCAGCGCGGCGAGCTGTGCGTCCGTGGCCCCATGATCATGAGCGGCTACGTTAACAACCCCGAGGCTACAAACGCTCTCATCGACAAGGACGGCTGGCTGCACAGCGGCGACATCGCCTACTGGGACGAGGACGAGCACTTCTTCATCGTGGACCGGCTGAAGAGCCTGATCAAATACAAGGGCTACCAGGTAGCCCCAGCCGAACTGGAGAGCATCCTGCTGCAACACCCCAACATCTTCGACGCCGGGGTCGCCGGCCTGCCCGACGACGATGCCGGCGAGCTGCCCGCCGCAGTCGTCGTGCTGGAACACGGTAAAACCATGACCGAGAAGGAGATCGTGGACTATGTGGCCAGCCAGGTTACAACCGCCAAGAAGCTGCGCGGTGGTGTTGTGTTCGTGGACGAGGTGCCTAAAGGACTGACCGGCAAGTTGGACGCCCGCAAGATCCGCGAGATTCTCATTAAGGCCAAGAAGGGCGGCAAGATCGCCGTGTAATAATCGTTCCTCTAGAGACGCGCAGGAGAAATTAATCAAGACTAGTACACTCCCCGTCGATCAGGGTGGTTACGTCAGTCACCGGTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAAT'
#print(consSeq)
alignments = pairwise2.align.globalxx(cons1, consSeq)
print(pairwise2.format_alignment(*alignments[0]))
'''
