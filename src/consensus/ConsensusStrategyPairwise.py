from consensus.ConsensusStrategy import ConsensusStrategy
from consensus.ReferenceConsensusGenerator import ReferenceConsensusGenerator
from Bio.Align import PairwiseAligner
from statistics import mean
from consensus.consensusStrategyPairwiseFunctions import identify_differences_from_indices, find_in_string_indices_of_character, inject_difference_into_sequence
from collections import Counter
import numpy as np

class ConsensusStrategyPairwise(ConsensusStrategy):
    def __init__(self):
        aligner = PairwiseAligner()
        aligner.mismatch_score = -1
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.5
        self.aligner = aligner

    def find_pairwise_score_and_all_differences_between_two_sequences(self, originalSequence, differentSequence):
        alignments = self.aligner.align(originalSequence, differentSequence)
        alignment = alignments[0]
        originalSequenceAlignment, indelIndicator, differentSequenceAlignment, _ = alignment.format().split("\n")
        differencesFromOriginal = []
        insertionIndices = find_in_string_indices_of_character(originalSequenceAlignment, "-")
        differencesFromOriginal.extend(identify_differences_from_indices("insertion", insertionIndices, originalSequenceAlignment, differentSequenceAlignment))

        deletionIndices = find_in_string_indices_of_character(differentSequenceAlignment, "-")
        differencesFromOriginal.extend(identify_differences_from_indices("deletion", deletionIndices, originalSequenceAlignment, ""))

        mutationIndices = find_in_string_indices_of_character(indelIndicator, ".")
        differencesFromOriginal.extend(identify_differences_from_indices("mutation", mutationIndices, originalSequenceAlignment, differentSequenceAlignment))

        return alignment.score, differencesFromOriginal

    def find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences(self, candidateSequence, readSequences):
        alignedScores = []
        allReadSequenceDifferences = []
        for readSequence in readSequences:
            score, readSequenceDifferences = self.find_pairwise_score_and_all_differences_between_two_sequences(candidateSequence, readSequence)
            alignedScores.append(score)
            allReadSequenceDifferences.extend(readSequenceDifferences)
        return mean(alignedScores), allReadSequenceDifferences

    def generate_consensus_sequence_from_biopython_records(self, binRecords: list) -> str:
        binSequences = [str(record.seq) for record in binRecords]
        referenceConsensusGenerator = ReferenceConsensusGenerator()
        referenceSequence = referenceConsensusGenerator.generate_consensus_sequence(binSequences)
        bestScore = -np.inf
        candidateSequence = referenceSequence[:]
        currentScore, currentDifferences = self.find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences(candidateSequence, binSequences)
        while currentScore > bestScore:
            if len(currentDifferences) == 0: return candidateSequence
            mostCommonDifference = Counter(currentDifferences).most_common(1)[0][0]
            nextSequence = inject_difference_into_sequence(candidateSequence, mostCommonDifference)
            bestScore = currentScore
            currentScore, currentDifferences = self.find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences(nextSequence, binSequences)
            if currentScore > bestScore: candidateSequence = nextSequence[:]
        return candidateSequence
