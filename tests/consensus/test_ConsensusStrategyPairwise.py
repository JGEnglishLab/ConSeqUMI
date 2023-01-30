from pytestConsensusFixtures import consensusSequence, readSequences, readSequenceRecords, simpleInsert, simpleString, middleInsertIndex, readSequenceDifferences
import pytest
import random
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)

from consensus import ConsensusStrategyPairwise
from consensus.ReferenceConsensusGenerator import ReferenceConsensusGenerator


def test__consensus_strategy_pairwise__initialization():
    pairwise = ConsensusStrategyPairwise.ConsensusStrategyPairwise()
    assert pairwise.aligner.mismatch_score == -1
    assert pairwise.aligner.open_gap_score == -1
    assert pairwise.aligner.extend_gap_score == -0.5

@pytest.fixture
def consensusStrategyPairwise():
    return ConsensusStrategyPairwise.ConsensusStrategyPairwise()

@pytest.fixture
def stretchLength():
    return 4

@pytest.fixture
def simpleStringWithThreeInserts(simpleInsert, simpleString, middleInsertIndex, stretchLength):
    simpleStringWithThreeInserts = simpleInsert * stretchLength + simpleString[:middleInsertIndex] + simpleInsert + simpleString[middleInsertIndex:] + simpleInsert * stretchLength
    return simpleStringWithThreeInserts

def test__consensus_strategy_pairwise__find_pairwise_score_and_all_differences_between_two_sequences__finds_all_inserts(simpleInsert, simpleString, middleInsertIndex, stretchLength, simpleStringWithThreeInserts, consensusStrategyPairwise):
    insertionDifferences = [
        (0,0,simpleInsert*stretchLength),
        (middleInsertIndex, middleInsertIndex, simpleInsert),
        (len(simpleString), len(simpleString), simpleInsert*stretchLength),
    ]
    _, insertionDifferencesOutput = consensusStrategyPairwise.find_pairwise_score_and_all_differences_between_two_sequences(simpleString, simpleStringWithThreeInserts)
    assert insertionDifferencesOutput == insertionDifferences

def test__consensus_strategy_pairwise__find_pairwise_score_and_all_differences_between_two_sequences__finds_all_deletions(simpleInsert, simpleString, middleInsertIndex, stretchLength, simpleStringWithThreeInserts, consensusStrategyPairwise):
    deletionDifferences = [
        (0,stretchLength,""),
        (middleInsertIndex+stretchLength, middleInsertIndex+stretchLength+len(simpleInsert), ""),
        (len(simpleStringWithThreeInserts)-stretchLength, len(simpleStringWithThreeInserts), ""),
    ]
    _, deletionDifferencesOutput = consensusStrategyPairwise.find_pairwise_score_and_all_differences_between_two_sequences(simpleStringWithThreeInserts, simpleString)
    assert deletionDifferencesOutput == deletionDifferences

def test__consensus_strategy_pairwise__find_pairwise_score_and_all_differences_between_two_sequences__finds_all_mutations(simpleInsert, simpleString, middleInsertIndex, stretchLength, simpleStringWithThreeInserts, consensusStrategyPairwise):
    simpleStringWithThreeMutations = simpleInsert + simpleString[1:middleInsertIndex] + simpleInsert + simpleString[middleInsertIndex+1:-1] + simpleInsert
    mutationDifferences = [
        (0,1,simpleInsert),
        (middleInsertIndex, middleInsertIndex+len(simpleInsert), simpleInsert),
        (len(simpleString)-1, len(simpleString), simpleInsert),
    ]
    _, mutationDifferencesOutput = consensusStrategyPairwise.find_pairwise_score_and_all_differences_between_two_sequences(simpleString, simpleStringWithThreeMutations)
    assert mutationDifferencesOutput == mutationDifferences

def test__consensus_strategy_pairwise__find_pairwise_score_and_all_differences_between_two_sequences__finds_all_deletions_after_insertion(simpleInsert, simpleString, middleInsertIndex, stretchLength, simpleStringWithThreeInserts, consensusStrategyPairwise):
    simpleStringWithLongerFrontInsert = simpleInsert*(stretchLength+1) + simpleString
    deletionDifferences = [
        (0,0,simpleInsert),
        (middleInsertIndex+stretchLength, middleInsertIndex+stretchLength+len(simpleInsert), ""),
        (len(simpleStringWithThreeInserts)-stretchLength, len(simpleStringWithThreeInserts), ""),
    ]
    _, deletionDifferencesOutput = consensusStrategyPairwise.find_pairwise_score_and_all_differences_between_two_sequences(simpleStringWithThreeInserts, simpleStringWithLongerFrontInsert)
    assert deletionDifferencesOutput == deletionDifferences

def calculate_average_pairwise_alignment_score_for_tests(consensusSequence):
    baseScore = len(consensusSequence)
    initialErrorScore = baseScore - 1
    allScores = []
    insertionScores = [baseScore-i*0.5-1 for i in range(5)]
    deletionScores = [baseScore-i*1.5-2 for i in range(5)]
    mutationScores = [baseScore - i*2 for i in range(1,5)]
    allScores.extend(insertionScores)
    allScores.extend(deletionScores)
    allScores.extend(mutationScores)
    return sum(allScores) / len(allScores)

def test__consensus_strategy_pairwise__find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences__finds_average_score(consensusSequence, readSequences, consensusStrategyPairwise):
    averagePairwiseAlignmentScore = calculate_average_pairwise_alignment_score_for_tests(consensusSequence)
    averagePairwiseAlignmentScoreOutput, _ = consensusStrategyPairwise.find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences(consensusSequence, readSequences)
    assert averagePairwiseAlignmentScore == averagePairwiseAlignmentScoreOutput

def test__consensus_strategy_pairwise__find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences__finds_all_differences(consensusSequence, readSequences, consensusStrategyPairwise, readSequenceDifferences):
    _,  readSequenceDifferencesOutput = consensusStrategyPairwise.find_average_pairwise_alignment_score_and_all_differences_between_candidate_sequence_and_binned_sequences(consensusSequence, readSequences)
    assert readSequenceDifferencesOutput == readSequenceDifferences

def test__consensus_strategy_pairwise__generate_consensus_sequence_from_biopython_records(consensusSequence, readSequences, consensusStrategyPairwise, readSequenceRecords):
    referenceConsensusGenerator = ReferenceConsensusGenerator()
    referenceSequence = referenceConsensusGenerator.generate_consensus_sequence(readSequences)
    consensusSequenceOutput = consensusStrategyPairwise.generate_consensus_sequence_from_biopython_records(readSequenceRecords)
    assert consensusSequenceOutput != referenceSequence
    assert consensusSequenceOutput == consensusSequence
