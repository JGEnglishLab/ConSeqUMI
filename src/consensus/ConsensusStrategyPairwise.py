from consensus.ConsensusStrategy import ConsensusStrategy
from Bio.Align import PairwiseAligner
from statistics import mean

class ConsensusStrategyPairwise(ConsensusStrategy):
    def __init__(self):
        aligner = PairwiseAligner()
        aligner.mismatch_score = -1
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.5
        self.aligner = aligner

    def find_average_pairwise_alignment_score(self, candidateSequence, readSequences):
        alignedScores = []
        for readSequence in readSequences:
            alignments = self.aligner.align(candidateSequence, readSequence)
            alignedScores.append(alignments[0].score)
        return mean(alignedScores)

    def find_all_differences_between_two_sequences(self, originalSequence, differentSequence):
        alignments = self.aligner.align(originalSequence, differentSequence)
        matchingIndices = alignments[0].aligned

        originalSequenceMatchingIndices, differentSequenceMatchingIndices = list(matchingIndices[0]), list(matchingIndices[1])

        startingIndices = (0,0)
        originalSequenceMatchingIndices.insert(0,startingIndices)
        differentSequenceMatchingIndices.insert(0,startingIndices)

        originalSequenceEndingIndices = (len(originalSequence),len(originalSequence))
        differentSequenceEndingIndices = (len(differentSequence),len(differentSequence))
        originalSequenceMatchingIndices.append(originalSequenceEndingIndices)
        differentSequenceMatchingIndices.append(differentSequenceEndingIndices)

        differencesFromOriginal = []
        numMatchingIndices = len(originalSequenceMatchingIndices)
        print()
        print(alignments[0].format().split("\n")[1])
        for i in range(numMatchingIndices - 1):
            originalSequenceDifferenceStartIndex = originalSequenceMatchingIndices[i][1]
            originalSequenceDifferenceEndIndex = originalSequenceMatchingIndices[i + 1][0]
            differentSequenceInsertStartIndex = differentSequenceMatchingIndices[i][1]
            differentSequenceInsertEndIndex = differentSequenceMatchingIndices[i + 1][0]
            differentSequenceInsert = differentSequence[differentSequenceInsertStartIndex:differentSequenceInsertEndIndex]
            if originalSequenceDifferenceStartIndex == originalSequenceDifferenceEndIndex and len(differentSequenceInsert) == 0: continue
            differencesFromOriginal.append((originalSequenceDifferenceStartIndex, originalSequenceDifferenceEndIndex, differentSequenceInsert))

        print(originalSequence)
        print(differentSequence)
        print(alignments[0])
        print(originalSequenceMatchingIndices)
        print(differentSequenceMatchingIndices)
        print(differencesFromOriginal)

        return differencesFromOriginal


    def generate_consensus_sequence_from_records(self, binRecords: list) -> str:
        pass
