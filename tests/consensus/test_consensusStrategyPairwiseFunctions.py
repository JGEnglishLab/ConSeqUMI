import pytest
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
from consensus import consensusStrategyPairwiseFunctions

@pytest.fixture
def simpleInsert():
    return "-"

@pytest.fixture
def simpleString():
    return "A"*10

@pytest.fixture
def simpleStringWithSingleFrontInsert(simpleInsert, simpleString):
    return simpleInsert + simpleString

@pytest.fixture
def singleFrontInsertIndices():
    return [[0]]

@pytest.fixture
def simpleStringWithSingleBackInsert(simpleInsert, simpleString):
    return simpleString + simpleInsert

@pytest.fixture
def singleBackInsertIndices(simpleStringWithSingleBackInsert):
    return [[len(simpleStringWithSingleBackInsert)-1]]

@pytest.fixture
def middleInsertIndex(simpleString):
    return len(simpleString) // 2

@pytest.fixture
def simpleStringWithSingleMiddleInsert(simpleInsert, simpleString, middleInsertIndex):
    return simpleString[:middleInsertIndex] + simpleInsert + simpleString[middleInsertIndex:]

@pytest.fixture
def singleMiddleInsertIndices(middleInsertIndex):
    return [[middleInsertIndex]]

@pytest.fixture
def singleBackMutationIndices(simpleString):
    return [[len(simpleString)-1]]

@pytest.fixture
def simpleStringWithSingleFrontMutation(simpleInsert, simpleString):
    return simpleInsert + simpleString[1:]

@pytest.fixture
def simpleStringWithSingleBackMutation(simpleInsert, simpleString):
    return simpleString[:-1] + simpleInsert

@pytest.fixture
def simpleStringWithSingleMiddleMutation(simpleInsert, simpleString, middleInsertIndex):
    return simpleString[:middleInsertIndex] + simpleInsert + simpleString[middleInsertIndex+1:]

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__for_single_instance_of_character_in_front(simpleInsert, simpleStringWithSingleFrontInsert, singleFrontInsertIndices):
    singleFrontInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=simpleStringWithSingleFrontInsert, character=simpleInsert)
    assert singleFrontInsertIndicesOutput == singleFrontInsertIndices

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__for_single_instance_of_character_in_back(simpleInsert, simpleStringWithSingleBackInsert, singleBackInsertIndices):
    singleBackInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=simpleStringWithSingleBackInsert, character=simpleInsert)
    assert singleBackInsertIndicesOutput == singleBackInsertIndices

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__for_single_instance_of_character_in_middle(simpleInsert, simpleStringWithSingleMiddleInsert, singleMiddleInsertIndices):
    singleMiddleInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=simpleStringWithSingleMiddleInsert, character=simpleInsert)
    assert singleMiddleInsertIndicesOutput == singleMiddleInsertIndices

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__for_multiple_instances_of_character_in_front(simpleInsert, simpleStringWithSingleFrontInsert):
    stringWithDoubleFrontInsert = simpleInsert + simpleStringWithSingleFrontInsert
    doubleFrontInsertIndices = [[0,1]]
    doubleFrontInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithDoubleFrontInsert, character=simpleInsert)
    assert doubleFrontInsertIndicesOutput == doubleFrontInsertIndices

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__for_multiple_instances_of_character_in_back(simpleInsert, simpleStringWithSingleBackInsert):
    stringWithDoubleBackInsert = simpleStringWithSingleBackInsert + simpleInsert
    doubleBackInsertIndices = [[len(stringWithDoubleBackInsert)-2, len(stringWithDoubleBackInsert)-1]]
    doubleBackInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithDoubleBackInsert, character=simpleInsert)
    assert doubleBackInsertIndicesOutput == doubleBackInsertIndices

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__for_multiple_instances_of_character_in_middle(simpleInsert, simpleStringWithSingleMiddleInsert, middleInsertIndex):
    stringWithDoubleMiddleInsert = simpleStringWithSingleMiddleInsert[:middleInsertIndex] + simpleInsert + simpleStringWithSingleMiddleInsert[middleInsertIndex:]
    middleInsertIndices = [[middleInsertIndex, middleInsertIndex+1]]
    middleInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithDoubleMiddleInsert, character=simpleInsert)
    assert middleInsertIndicesOutput == middleInsertIndices

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__with_separate_single_and_multiple_instances_of_character(simpleInsert, simpleString):
    stringWithDoubleFrontInsert = simpleInsert * 2 + simpleString
    stringWithDoubleFrontInsertAndSingleBackInsert = stringWithDoubleFrontInsert + simpleInsert
    indices = [[0,1], [len(stringWithDoubleFrontInsertAndSingleBackInsert)-1]]
    indicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithDoubleFrontInsertAndSingleBackInsert, character=simpleInsert)
    assert indicesOutput == indices

@pytest.fixture
def singleFrontInsertDifference(simpleInsert):
    return (0, 0, simpleInsert)

@pytest.fixture
def singleBackInsertDifference(simpleInsert, simpleString):
    return (len(simpleInsert + simpleString)-1, len(simpleInsert + simpleString)-1, simpleInsert)

@pytest.fixture
def singleMiddleInsertDifference(simpleInsert, middleInsertIndex):
    return (middleInsertIndex, middleInsertIndex, simpleInsert)

@pytest.fixture
def singleFrontDeletionDifference():
    return (0, 1, "")

@pytest.fixture
def singleBackDeletionDifference(simpleStringWithSingleBackInsert):
    return (len(simpleStringWithSingleBackInsert)-1,len(simpleStringWithSingleBackInsert),"")

@pytest.fixture
def singleMiddleDeletionDifference(middleInsertIndex):
    return (middleInsertIndex,middleInsertIndex+1,"")

@pytest.fixture
def singleFrontMutationDifference(simpleInsert):
    return (0, 1, simpleInsert)

@pytest.fixture
def singleBackMutationDifference(simpleInsert, simpleString):
    return (len(simpleString)-1,len(simpleString), simpleInsert)

@pytest.fixture
def singleMiddleMutationDifference(simpleInsert, middleInsertIndex):
    return (middleInsertIndex,middleInsertIndex+1,simpleInsert)

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_insertion_to_front(simpleString, simpleStringWithSingleFrontInsert, singleFrontInsertDifference):
    simpleStringWithSingleFrontInsertOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, singleFrontInsertDifference)
    assert simpleStringWithSingleFrontInsertOutput == simpleStringWithSingleFrontInsert

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_insertion_to_back(simpleString, simpleStringWithSingleBackInsert, singleBackInsertDifference):
    simpleStringWithSingleBackInsertOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, singleBackInsertDifference)
    assert simpleStringWithSingleBackInsertOutput == simpleStringWithSingleBackInsert

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_insertion_to_middle(simpleString, simpleStringWithSingleMiddleInsert, singleMiddleInsertDifference):
    simpleStringWithSingleMiddleInsertOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, singleMiddleInsertDifference)
    assert simpleStringWithSingleMiddleInsertOutput == simpleStringWithSingleMiddleInsert

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_deletion_to_front(simpleString, simpleStringWithSingleFrontInsert, singleFrontDeletionDifference):
    simpleStringOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleStringWithSingleFrontInsert, singleFrontDeletionDifference)
    assert simpleStringOutput == simpleString

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_deletion_to_back(simpleString, simpleStringWithSingleBackInsert, singleBackDeletionDifference):
    simpleStringOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleStringWithSingleBackInsert, singleBackDeletionDifference)
    assert simpleStringOutput == simpleString

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_deletion_to_middle(simpleString, simpleStringWithSingleMiddleInsert, singleMiddleDeletionDifference):
    simpleStringOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleStringWithSingleMiddleInsert, singleMiddleDeletionDifference)
    assert simpleStringOutput == simpleString

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_mutation_to_front(simpleInsert, simpleString, singleFrontMutationDifference, simpleStringWithSingleFrontMutation):
    simpleStringWithSingleFrontMutationOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, singleFrontMutationDifference)
    assert simpleStringWithSingleFrontMutationOutput == simpleStringWithSingleFrontMutation

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_mutation_to_back(simpleInsert, simpleString, singleBackMutationDifference, simpleStringWithSingleBackMutation):
    simpleStringWithSingleBackMutationOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, singleBackMutationDifference)
    assert simpleStringWithSingleBackMutationOutput == simpleStringWithSingleBackMutation

def test__consensus_strategy_pairwise_functions__inject_difference_into_sequence__inject_mutation_to_middle(simpleInsert, simpleString, middleInsertIndex, singleMiddleMutationDifference, simpleStringWithSingleMiddleMutation):
    simpleStringWithSingleMiddleMutationOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, singleMiddleMutationDifference)
    assert simpleStringWithSingleMiddleMutationOutput == simpleStringWithSingleMiddleMutation

@pytest.fixture
def simpleStringWithFrontAndBackInsert(simpleStringWithSingleFrontInsert, simpleInsert):
    return simpleStringWithSingleFrontInsert + simpleInsert

@pytest.fixture
def frontAndBackIndices(singleFrontInsertIndices, simpleStringWithFrontAndBackInsert):
    return singleFrontInsertIndices + [[len(simpleStringWithFrontAndBackInsert)-1]]

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_insertion_in_front(simpleString, singleFrontInsertIndices, simpleStringWithSingleFrontInsert, singleFrontInsertDifference):
    singleFrontInsertDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("insertion", singleFrontInsertIndices, simpleStringWithSingleFrontInsert)
    assert singleFrontInsertDifferenceOutput == [singleFrontInsertDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, *singleFrontInsertDifferenceOutput) == simpleStringWithSingleFrontInsert

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_insertion_indices_in_back(simpleString, singleBackInsertIndices, simpleStringWithSingleBackInsert, singleBackInsertDifference):
    singleBackInsertDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("insertion", singleBackInsertIndices, simpleStringWithSingleBackInsert)
    assert singleBackInsertDifferenceOutput == [singleBackInsertDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, *singleBackInsertDifferenceOutput) == simpleStringWithSingleBackInsert

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_insertion_indices_in_middle(simpleString, singleMiddleInsertIndices, simpleStringWithSingleMiddleInsert, singleMiddleInsertDifference):
    singleMiddleInsertDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("insertion", singleMiddleInsertIndices, simpleStringWithSingleMiddleInsert)
    assert singleMiddleInsertDifferenceOutput == [singleMiddleInsertDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, *singleMiddleInsertDifferenceOutput) == simpleStringWithSingleMiddleInsert

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_insertion_indices_in_front_and_back(simpleInsert, simpleStringWithSingleFrontInsert, singleFrontInsertDifference, simpleStringWithFrontAndBackInsert, frontAndBackIndices):
    frontAndBackDifferences = [singleFrontInsertDifference, (len(simpleStringWithFrontAndBackInsert)-1, len(simpleStringWithFrontAndBackInsert)-1, simpleInsert)]
    frontAndBackDifferencesOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("insertion", frontAndBackIndices, simpleStringWithFrontAndBackInsert)
    assert frontAndBackDifferencesOutput == frontAndBackDifferences

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_deletion_indices_in_front(simpleString, singleFrontInsertIndices, singleFrontDeletionDifference, simpleStringWithSingleFrontInsert):
    singleFrontDeletionDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("deletion", singleFrontInsertIndices)
    assert singleFrontDeletionDifferenceOutput == [singleFrontDeletionDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleStringWithSingleFrontInsert, *singleFrontDeletionDifferenceOutput) == simpleString

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_deletion_indices_in_back(simpleString, singleBackInsertIndices, singleBackDeletionDifference, simpleStringWithSingleBackInsert):
    singleBackDeletionDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("deletion", singleBackInsertIndices)
    assert singleBackDeletionDifferenceOutput == [singleBackDeletionDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleStringWithSingleBackInsert, *singleBackDeletionDifferenceOutput) == simpleString

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_deletion_indices_in_middle(simpleString, singleMiddleInsertIndices, singleMiddleDeletionDifference, simpleStringWithSingleMiddleInsert):
    singleMiddleDeletionDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("deletion", singleMiddleInsertIndices)
    assert singleMiddleDeletionDifferenceOutput == [singleMiddleDeletionDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleStringWithSingleMiddleInsert, *singleMiddleDeletionDifferenceOutput) == simpleString

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_deletion_indices_in_front_and_back(simpleInsert, simpleStringWithSingleFrontInsert, singleFrontDeletionDifference, simpleStringWithFrontAndBackInsert, frontAndBackIndices):
    frontAndBackDifferences = [singleFrontDeletionDifference, (len(simpleStringWithFrontAndBackInsert)-1, len(simpleStringWithFrontAndBackInsert), "")]
    frontAndBackDifferencesOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("deletion", frontAndBackIndices)
    assert frontAndBackDifferencesOutput == frontAndBackDifferences

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_mutation_indices_in_front(simpleString, singleFrontInsertIndices, simpleStringWithSingleFrontMutation, singleFrontMutationDifference):
    singleFrontMutationDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("mutation", singleFrontInsertIndices, simpleStringWithSingleFrontMutation)
    assert singleFrontMutationDifferenceOutput == [singleFrontMutationDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, *singleFrontMutationDifferenceOutput) == simpleStringWithSingleFrontMutation

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_mutation_indices_in_back(simpleString, singleBackMutationIndices, simpleStringWithSingleBackMutation, singleBackMutationDifference):
    singleBackMutationDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("mutation", singleBackMutationIndices, simpleStringWithSingleBackMutation)
    assert singleBackMutationDifferenceOutput == [singleBackMutationDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, *singleBackMutationDifferenceOutput) == simpleStringWithSingleBackMutation

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_mutation_indices_in_middle(simpleString, singleMiddleInsertIndices, simpleStringWithSingleMiddleMutation, singleMiddleMutationDifference):
    singleMiddleMutationDifferenceOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("mutation", singleMiddleInsertIndices, simpleStringWithSingleMiddleMutation)
    assert singleMiddleMutationDifferenceOutput == [singleMiddleMutationDifference]
    assert consensusStrategyPairwiseFunctions.inject_difference_into_sequence(simpleString, *singleMiddleMutationDifferenceOutput) == simpleStringWithSingleMiddleMutation

def test__consensus_strategy_pairwise_functions__identify_differences_from_indices__for_mutation_indices_in_front_and_back(simpleInsert, simpleStringWithSingleFrontInsert, singleFrontMutationDifference, simpleStringWithFrontAndBackInsert, frontAndBackIndices):
    frontAndBackDifferences = [singleFrontMutationDifference, (len(simpleStringWithFrontAndBackInsert)-1, len(simpleStringWithFrontAndBackInsert), simpleInsert)]
    frontAndBackDifferencesOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("mutation", frontAndBackIndices, simpleStringWithFrontAndBackInsert)
    assert frontAndBackDifferencesOutput == frontAndBackDifferences
