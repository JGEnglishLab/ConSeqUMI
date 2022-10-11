import numpy as np
import pandas as pd
import re
from collections import defaultdict, Counter
from statistics import median, mean
from Bio import SeqIO
from Bio.Align import PairwiseAligner

def find_average_aligned_score(candidateSeq, binSeqs):
    aligner = PairwiseAligner()
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    scores = []
    for binSeq in binSeqs:
        alignments = aligner.align(candidateSeq, binSeq)
        scores.append(alignments[0].score)
    return mean(scores)

def find_aligned_differences(seq1, seq2):

    aligner = PairwiseAligner()
    alignments = aligner.align(seq1,seq2)
    indices = alignments[0].aligned
    a_idx, b_idx = list(indices[0]), list(indices[1])
    a_idx.insert(0, (0,0))
    b_idx.insert(0, (0,0))
    a_idx.append((len(seq1),len(seq1)))
    b_idx.append((len(seq2),len(seq2)))

    diffs = []

    for i in range(len(a_idx)-1):
        start = a_idx[i][1]
        end = a_idx[i+1][0]
        insert = seq2[b_idx[i][1]:b_idx[i+1][0]]
        if start == end and len(insert)==0: continue
        diffs.append((start, end, insert))
    return diffs

def update_candidate_seq_by_common_diffs(candidateSeq, binSeqs, cutoff_percent = None):
    diffs = []
    for binSeq in binSeqs: diffs.extend(find_aligned_differences(candidateSeq, binSeq))
    finalSeq = candidateSeq[:]

    if cutoff_percent:
        cutoff = len(binSeqs) * cutoff_percent
        commonDiffs = sorted([key for key,value in Counter(diffs).items() if value > cutoff],reverse=True)
        for diff in commonDiffs:
            start, end, insert = diff
            finalSeq = finalSeq[:start] + insert + finalSeq[end:]
    else:
        most_common_diff = Counter(diffs).most_common(1)[0][0]
        start, end, insert = most_common_diff
        finalSeq = finalSeq[:start] + insert + finalSeq[end:]
    return finalSeq


class ConsensusSequenceGenerator:
    def __init__(self, match_reward=1, mismatch_penalty=-1, gap_penalty=-2):
        self.GAP_PENALTY = gap_penalty
        self.SAME_AWARD = match_reward
        self.DIFFERENCE_PENALTY = mismatch_penalty

    def generate_consensus_sequence_from_file(self, binFile, cutoff_percent=None):
        binSeqs = [str(record.seq) for record in SeqIO.parse(binFile, "fastq")]
        return self.generate_consensus_sequence(binSeqs, cutoff_percent)

    def generate_consensus_sequence(self, binSeqs, cutoff_percent):
        diffs = []
        refSeq = find_reference_sequence(binSeqs)
        if cutoff_percent: return update_candidate_seq_by_common_diffs(refSeq, binSeqs, cutoff_percent)
        bestScore = -np.inf
        candidateSeq = refSeq[:]
        curScore = find_average_aligned_score(candidateSeq, binSeqs)
        while curScore >= bestScore:
            bestScore = curScore
            tempSeq = update_candidate_seq_by_common_diffs(candidateSeq, binSeqs)
            curScore = find_average_aligned_score(tempSeq, binSeqs)
            if curScore >= bestScore: candidateSeq = tempSeq
        return candidateSeq

        '''
        cutoff_percents = [round(i,2) for i in np.arange(0.3, 0.61, 0.01)]
        seqs = []
        for c in cutoff_percents:
            cutoff2 = len(binSeqs) * c
            commonDiffs2 = sorted([key for key,value in Counter(diffs).items() if value > cutoff2],reverse=True)
            finalSeq2 = refSeq[:]
            for diff in commonDiffs2:
                start, end, insert = diff
                finalSeq2 = finalSeq2[:start] + insert + finalSeq2[end:]
            seqs.append(finalSeq2)
        print(",".join(seqs))
        '''
        return finalSeq, refSeq

    def __find_alignment_matrices__(self, seq_a, seq_b):

        len_a = len(seq_a)
        len_b = len(seq_b)
        alignment_matrix = np.empty((len_a + 1, len_b + 1))
        arrow_matrix = np.zeros((len_a + 1, len_b + 1, 3), dtype=bool)

        alignment_matrix[0, 0] = 0

        for i in range(1, len_a + 1):
            alignment_matrix[i, 0] = alignment_matrix[i - 1, 0]\
                                     + self.GAP_PENALTY
            arrow_matrix[i, 0, 1] = True

        for j in range(1, len_b + 1):
            alignment_matrix[0, j] = alignment_matrix[0, j - 1]\
                                     + self.GAP_PENALTY
            arrow_matrix[0, j, 0] = True

        for i in range(1, len_a + 1):
            for j in range(1, len_b + 1):

                right_score = alignment_matrix[i, j - 1] + self.GAP_PENALTY
                down_score = alignment_matrix[i - 1, j] + self.GAP_PENALTY
                if seq_a[i - 1] == seq_b[j - 1]:
                    diag_score = alignment_matrix[i - 1, j - 1]\
                                 + self.SAME_AWARD
                else:
                    diag_score = alignment_matrix[i - 1, j - 1] \
                                 + self.DIFFERENCE_PENALTY

                best_score = np.max((right_score, down_score, diag_score))

                arrow_matrix[i, j, 0] = right_score == best_score
                arrow_matrix[i, j, 1] = down_score == best_score
                arrow_matrix[i, j, 2] = diag_score == best_score

                alignment_matrix[i, j] = best_score

        return alignment_matrix, arrow_matrix

    def __find_string_differences__(self, seq_b, arrow_matrix):
        i, j, _ = arrow_matrix.shape
        i -= 1
        j -= 1
        string_differences = ['' for _ in range(i*2+1)]
        while i > 0 or j > 0:
            right_move, down_move, diag_move = arrow_matrix[i,j]
            seq_b_val = seq_b[j-1]
            if right_move:
                altered_index = i*2
                string_differences[altered_index] = seq_b_val + string_differences[altered_index]
                j -=1
                continue
            if down_move: i -= 1; continue
            if diag_move:
                altered_index = (i-1)*2+1
                string_differences[altered_index] += seq_b_val
                i -= 1
                j -= 1
                continue
        return string_differences

def adjust_all_string_lengths(strs, buffer_length):
    max_length = len(max(strs, key = len))
    strs = [x.ljust(max_length+buffer_length) for x in strs]
    return strs

def initialize_consensus_output_string(strs, excerpt_length):
    tempList = [x[:excerpt_length] for x in strs]
    final = max(set(tempList), key = tempList.count)
    return final

def find_next_character(subStrings, tempPattern):
    pattern = re.compile(tempPattern + '.')
    baseCountDict = defaultdict(list)
    all_characters = []
    for x in subStrings:
        r = pattern.search(x)
        if r: baseCountDict[r.group(0)[-1]].append(r.start(0)); all_characters.append(r.group(0)[-1])
    c = Counter(all_characters)
    percent_benchmark = .2
    benchmark = int(len(subStrings)*percent_benchmark)
    common = c.most_common()
    if len(common) > 1 and common[1][1] > benchmark:
        min_val = min([median(baseCountDict[x[0]]) for x in common[:2]])
        res = [x[0] for x in common[:2] if median(baseCountDict[x[0]]) == min_val]
        return res[0]
    return c.most_common()[0][0]

def find_reference_sequence(strs, excerpt_length = 10, buffer_length = 20):
    strs = adjust_all_string_lengths(strs, buffer_length)
    final = initialize_consensus_output_string(strs, excerpt_length)
    for i in range(len(strs[0])-buffer_length):
        subStrings = [x[i:i+buffer_length] for x in strs]
        tempPattern = final[-excerpt_length:]
        nex = find_next_character(subStrings, tempPattern)
        final += nex
    return final.strip(' ')

if __name__ == "__main__":
    nw = ConsensusSequenceGenerator()
    seq1 =  'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'
    seqs = ['ATCGATCGATCGATCGATCGATCGAACGATCGATCGATCGATCGAACGATCCATCGATCGATCG',
            'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGTGATCGATCGATCG',
            'ATCGATCGATCGATCGATCGATCGATCGATCGATCCCCGATCGATCGATCGATCGATCGATCGATCG',
            'ATCGATCGATCGATCGATCGATCGATCGATTCGATCGATCGATCGATCGATCGATCGATCG',
            'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCC']
    consSeq = nw.generate_consensus_sequence(seqs)
    assert seq1 == consSeq
