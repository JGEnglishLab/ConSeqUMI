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
        self.binSize = 20
        self.consLength = 10

    def generate_consensus_sequence(self, file):

        seqStrs = [str(record.seq) for record in SeqIO.parse(file, "fastq")]
        self.numBins = round(median(map(len, seqStrs))/self.binSize)
        consWindows = []
        for seq in seqStrs:
            windows = [seq[i:i + self.consLength] for i in range(len(seq) - self.consLength + 1)]
            consWindows.append(np.array(windows))


        window = self.find_consensus_windows(seqStrs, consWindows, offset = 0)
        return window

    def return_possible_consensus_indices(self, l, offset):
        """Yield n number of sequential chunks from l."""
        d, r = divmod(len(l), self.numBins)
        numBinsRange = range(self.numBins)
        if offset: numBinsRange = range(self.numBins - 1)
        for i in numBinsRange:
            si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
            diff = (d+1 if i < r else d)
            yield np.array(range(si, si + diff - self.consLength + 1)) + offset


    def find_consensus_windows(self, seqStrs, seqWindows, offset = 0):

        binnedConsWindows = []
        allBinIndices = []
        remove = []
        for i in range(len(seqStrs)):
            seq = seqStrs[i]
            windows = seqWindows[i]
            if len(windows) / self.numBins < self.consLength: remove.append(i); continue
            tempIndices = [x for x in self.return_possible_consensus_indices(windows, offset)]
            tempSets = [set(windows[x]) for x in tempIndices]
            binnedConsWindows.append(tempSets)

            allBinIndices.append(tempIndices)
            #binSeq = [seq[x[0]:x[-1]+10] for x in indices]
            #binSeqs.append(binSeq)
        for i in remove[::-1]: del seqStrs[i]
        consSeqs = []
        #print(binSeqs[0][0])
        #print(binnedConsWindows[0][0])
        numBinsRange = range(self.numBins)
        if offset: numBinsRange = range(self.numBins - 1)
        for binIdx in numBinsRange:
            allPossibleCons = [list(x[binIdx]) for x in binnedConsWindows]
            allPossibleCons = list(itertools.chain.from_iterable(allPossibleCons))
            consSeqs.append(most_common(allPossibleCons))

        gapSeqPossibleCons = [[] for x in range(len(consSeqs)+1)]

        for i in range(len(seqStrs)):
            seq = seqStrs[i]
            binIndices = allBinIndices[i]
            newIndices = []
            if i == -1: print(seqStrs[i])
            for j in range(len(consSeqs)):
                idx = binIndices[j]
                binSeq = seq[idx[0]:idx[-1]+self.consLength]
                consSeq = consSeqs[j]
                alignments = pairwise2.align.localms(binSeq, consSeq, 2, -1, -1, -.5)
                alignments.sort(key=lambda x: x.start)
                if len(alignments) == 0: start = 0; end = self.consLength
                else: start = idx[0] + alignments[0].start; end = idx[0] + alignments[0].end
                if i == -1:
                    #print((idx[0],idx[-1]+self.consLength))
                    print(alignments)
                    print(binSeq)
                    print(consSeq)
                    print(binnedConsWindows[i][j])
                    print()
                newIndices.append((start, end))
                if j == 0: gapSeqPossibleCons[j].append(seq[:newIndices[j][0]])
                else: gapSeqPossibleCons[j].append(seq[newIndices[j-1][1]:newIndices[j][0]])
                if j == len(consSeqs)-1: gapSeqPossibleCons[j+1].append(seq[newIndices[j][1]:])

        gapSeqs = []
        for binIdx in numBinsRange:
            binAllPossibleGaps = gapSeqPossibleCons[binIdx]
            #binAllPossibleGaps = itertools.chain.from_iterable(binAllPossibleGaps)
            gapSeqs.append(most_common(binAllPossibleGaps))

        gapSeqs.append(most_common(gapSeqPossibleCons[-1]))

        allConsSeqs = []
        for i in range(len(consSeqs)):
            allConsSeqs.append(gapSeqs[i])
            allConsSeqs.append(consSeqs[i])
        allConsSeqs.append(gapSeqs[-1])
        #print('.'.join(allConsSeqs))
        output = ''.join(allConsSeqs)
        #print(output)
        while len(output) != 0 and output[-1] == 'A': output = output[:-1]
        return output

'''
binFile = 'scratchpaper_code/seq_bin194.fq'
csg = ConsensusSequenceGenerator()
consSeq = csg.generate_consensus_sequence(binFile)
cons1 = 'TGCCACCCGGCTTCAACGAGTACGACTTCGTGCCCGAGAGCTTCGACCGGGACAAAACCATCGCCCTGATCATGAACAGTAGTGGCAGTACCGGATTGCCCAAGGGCGTAGCCCTACCGCACCGCACCGCTTGTGTCCGATTCAGTCATGCCCGCGACCCCATCTTCGGCAACCAGATCATCCCCGACACCGCTATCCTCAGCGTGGTGCCATTTCACCACGGCTTCGGCATGTTCACCACGCTGGGCTACTTGATCTGCGGCTTTCGGGTCGTGCTCATGTACCGCTTCGAGGAGGAGCTATTCTTGCGCAGCTTGCAAGACTATAAGATTCAATCTGCCCTGCTGGTGCCCACACTATTTAGCTTCTTCGCTAAGAGCACTCTCATCGACAAGTACGACCTAAGCAACTTGCACGAGATCGCCAGCGGCGGGGCGCCGCTCAGCAAGGAGGTAGGTGAGGCCGTGGCCAAACGCTTCCACCTACCAGGCATCCGCCAGGGCTACGGCCTGACAGAAACAACCAGCGCCATTCTGATCACCCCCGAAGGGGACGACAAGCCTGGCGCAGTAGGCAAGGTGGTGCCCTTCTTCGAGGCTAAGGTGGTGGACTTGGACACCGGCAAGACACTGGGTGTGAACCAGCGCGGCGAGCTGTGCGTCCGTGGCCCCATGATCATGAGCGGCTACGTTAACAACCCCGAGGCTACAAACGCTCTCATCGACAAGGACGGCTGGCTGCACAGCGGCGACATCGCCTACTGGGACGAGGACGAGCACTTCTTCATCGTGGACCGGCTGAAGAGCCTGATCAAATACAAGGGCTACCAGGTAGCCCCAGCCGAACTGGAGAGCATCCTGCTGCAACACCCCAACATCTTCGACGCCGGGGTCGCCGGCCTGCCCGACGACGATGCCGGCGAGCTGCCCGCCGCAGTCGTCGTGCTGGAACACGGTAAAACCATGACCGAGAAGGAGATCGTGGACTATGTGGCCAGCCAGGTTACAACCGCCAAGAAGCTGCGCGGTGGTGTTGTGTTCGTGGACGAGGTGCCTAAAGGACTGACCGGCAAGTTGGACGCCCGCAAGATCCGCGAGATTCTCATTAAGGCCAAGAAGGGCGGCAAGATCGCCGTGTAATAATCGTTCCTCTAGAGACGCGCAGGAGAAATTAATCAAGACTAGTACACTCCCCGTCGATCAGGGTGGTTACGTCAGTCACCGGTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAAT'
#print(consSeq)
alignments = pairwise2.align.globalxx(cons1, consSeq)
print('\n')
print(pairwise2.format_alignment(*alignments[0]))
#'''
