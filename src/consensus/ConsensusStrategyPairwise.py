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

    def generate_consensus_sequence_from_records(self, binRecords: list) -> str:
        pass
