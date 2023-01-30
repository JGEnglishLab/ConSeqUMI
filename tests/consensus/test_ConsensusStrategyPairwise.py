from pytestConsensusFixtures import consensusSequence, readSequences, readSequenceRecords, simpleInsert, simpleString, middleInsertIndex
import pytest
import random
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)

from consensus import ConsensusStrategyPairwise


def test__consensus_strategy_pairwise__initialization():
    pairwise = ConsensusStrategyPairwise.ConsensusStrategyPairwise()
    assert pairwise.aligner.mismatch_score == -1
    assert pairwise.aligner.open_gap_score == -1
    assert pairwise.aligner.extend_gap_score == -0.5

@pytest.fixture
def consensusStrategyPairwise():
    return ConsensusStrategyPairwise.ConsensusStrategyPairwise()

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

def test__consensus_strategy_pairwise__find_average_pairwise_alignment_score(consensusSequence, readSequenceRecords, consensusStrategyPairwise):
    readSequences = [str(readSequenceRecord.seq) for readSequenceRecord in readSequenceRecords]
    averagePairwiseAlignmentScore = calculate_average_pairwise_alignment_score_for_tests(consensusSequence)
    averagePairwiseAlignmentScoreOutput = consensusStrategyPairwise.find_average_pairwise_alignment_score(consensusSequence, readSequences)
    assert averagePairwiseAlignmentScore == averagePairwiseAlignmentScoreOutput

@pytest.fixture
def stretchLength():
    return 4

@pytest.fixture
def simpleStringWithThreeInserts(simpleInsert, simpleString, middleInsertIndex, stretchLength):
    simpleStringWithThreeInserts = simpleInsert * stretchLength + simpleString[:middleInsertIndex] + simpleInsert + simpleString[middleInsertIndex:] + simpleInsert * stretchLength
    return simpleStringWithThreeInserts

def test__consensus_strategy_pairwise__find_all_differences_between_two_sequences__finds_all_inserts(simpleInsert, simpleString, middleInsertIndex, stretchLength, simpleStringWithThreeInserts, consensusStrategyPairwise):
    insertionDifferences = [
        (0,0,simpleInsert*stretchLength),
        (middleInsertIndex, middleInsertIndex, simpleInsert),
        (len(simpleString), len(simpleString), simpleInsert*stretchLength),
    ]
    insertionDifferencesOutput = consensusStrategyPairwise.find_all_differences_between_two_sequences(simpleString, simpleStringWithThreeInserts)
    assert insertionDifferencesOutput == insertionDifferences

def test__consensus_strategy_pairwise__find_all_differences_between_two_sequences__finds_all_deletions(simpleInsert, simpleString, middleInsertIndex, stretchLength, simpleStringWithThreeInserts, consensusStrategyPairwise):
    deletionDifferences = [
        (0,stretchLength,""),
        (middleInsertIndex+stretchLength, middleInsertIndex+stretchLength+len(simpleInsert), ""),
        (len(simpleStringWithThreeInserts)-stretchLength, len(simpleStringWithThreeInserts), ""),
    ]
    deletionDifferencesOutput = consensusStrategyPairwise.find_all_differences_between_two_sequences(simpleStringWithThreeInserts, simpleString)
    assert deletionDifferencesOutput == deletionDifferences

def test__consensus_strategy_pairwise__find_all_differences_between_two_sequences__finds_all_mutations(simpleInsert, simpleString, middleInsertIndex, stretchLength, simpleStringWithThreeInserts, consensusStrategyPairwise):
    simpleStringWithThreeMutations = simpleInsert + simpleString[1:middleInsertIndex] + simpleInsert + simpleString[middleInsertIndex+1:-1] + simpleInsert
    mutationDifferences = [
        (0,1,simpleInsert),
        (middleInsertIndex, middleInsertIndex+len(simpleInsert), simpleInsert),
        (len(simpleString)-1, len(simpleString), simpleInsert),
    ]
    mutationDifferencesOutput = consensusStrategyPairwise.find_all_differences_between_two_sequences(simpleString, simpleStringWithThreeMutations)
    assert mutationDifferencesOutput == mutationDifferences

def test__consensus_strategy_pairwise__find_all_differences_between_two_sequences__finds_all_deletions_after_insertion(simpleInsert, simpleString, middleInsertIndex, stretchLength, simpleStringWithThreeInserts, consensusStrategyPairwise):
    simpleStringWithLongerFrontInsert = simpleInsert*(stretchLength+1) + simpleString
    deletionDifferences = [
        (0,0,simpleInsert),
        (middleInsertIndex+stretchLength, middleInsertIndex+stretchLength+len(simpleInsert), ""),
        (len(simpleStringWithThreeInserts)-stretchLength, len(simpleStringWithThreeInserts), ""),
    ]
    deletionDifferencesOutput = consensusStrategyPairwise.find_all_differences_between_two_sequences(simpleStringWithThreeInserts, simpleStringWithLongerFrontInsert)
    assert deletionDifferencesOutput == deletionDifferences
