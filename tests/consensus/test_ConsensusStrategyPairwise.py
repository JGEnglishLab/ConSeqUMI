from pytestConsensusFixtures import consensusSequence, readSequences, readSequenceRecords
import pytest
import sys

sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from consensus.ConsensusStrategyPairwise import ConsensusStrategyPairwise

def test_consensus_strategy_pairwise_initialization():
    pairwise = ConsensusStrategyPairwise()
    assert pairwise.aligner.mismatch_score == -1
    assert pairwise.aligner.open_gap_score == -1
    assert pairwise.aligner.extend_gap_score == -0.5

@pytest.fixture
def consensusStrategyPairwise():
    return ConsensusStrategyPairwise()

def calculate_average_pairwise_alignment_score_for_tests(consensusSequence):
    baseScore = len(consensusSequence)
    initialErrorScore = baseScore - 1
    allScores = []
    insertionScores = [initialErrorScore-i*0.5 for i in range(5)]
    deletionScores = [initialErrorScore-i*1.5-1 for i in range(5)]
    mutationScores = [initialErrorScore-1 for i in range(5)]
    allScores.extend(insertionScores)
    allScores.extend(deletionScores)
    allScores.extend(mutationScores)
    return sum(allScores) / len(allScores)

def test_consensus_strategy_find_average_pairwise_alignment_score(consensusSequence, readSequenceRecords, consensusStrategyPairwise):
    readSequences = [str(readSequenceRecord.seq) for readSequenceRecord in readSequenceRecords]
    averagePairwiseAlignmentScore = calculate_average_pairwise_alignment_score_for_tests(consensusSequence)
    averagePairwiseAlignmentScoreOutput = consensusStrategyPairwise.find_average_pairwise_alignment_score(consensusSequence, readSequences)
    assert averagePairwiseAlignmentScore == averagePairwiseAlignmentScoreOutput

def test_consensus_strategy_generate_consensus_sequence_from_records(): pass
