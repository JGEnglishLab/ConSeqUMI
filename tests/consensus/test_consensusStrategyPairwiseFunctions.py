import pytest
import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
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

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__for_single_instance_of_character_in_front(simpleInsert, simpleStringWithSingleFrontInsert, singleFrontInsertIndices):
    singleFrontInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=simpleStringWithSingleFrontInsert, character=simpleInsert)
    assert singleFrontInsertIndicesOutput == singleFrontInsertIndices

@pytest.fixture
def simpleStringWithSingleBackInsert(simpleInsert, simpleString):
    return simpleString + simpleInsert

@pytest.fixture
def singleBackInsertIndices(simpleStringWithSingleBackInsert):
    return [[len(simpleStringWithSingleBackInsert)-1]]

def test__consensus_strategy_pairwise_functions__find_in_string_indices_of_character__for_single_instance_of_character_in_back(simpleInsert, simpleStringWithSingleBackInsert, singleBackInsertIndices):
    singleBackInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=simpleStringWithSingleBackInsert, character=simpleInsert)
    assert singleBackInsertIndicesOutput == singleBackInsertIndices

@pytest.fixture
def middleInsertIndex(simpleString):
    return len(simpleString) // 2

@pytest.fixture
def simpleStringWithSingleMiddleInsert(simpleInsert, simpleString, middleInsertIndex):
    return simpleString[:middleInsertIndex] + simpleInsert + simpleString[middleInsertIndex:]

@pytest.fixture
def singleMiddleInsertIndices(middleInsertIndex):
    return [[middleInsertIndex]]

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

# Continue from here

def test_consensus_strategy_pairwise_functions_format_difference_from_indices_for_insertion_indices():
    character = "A"
    fillerLength = 5
    alignment = character + "-" * fillerLength + character + "-" * fillerLength + character
    indices = [[0], [fillerLength+1], [len(alignment)-1]]
    differences = [
        (0, 0, character),
        (fillerLength+1, fillerLength+1, character),
        (len(alignment)-1, len(alignment)-1, character),
    ]
    differencesOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("insertion", indices, alignment)
    assert differencesOutput == differences

def test_consensus_strategy_pairwise_functions_format_differences_from_indices_for_deletion_indices():
    indices = [[0]]
    differences = [
        (0, 1, ""),
    ]
    differencesOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("deletion", indices)
    assert differencesOutput == differences

def test_consensus_strategy_pairwise_functions_format_difference_from_indices_for_mutation_indices():
    character = "A"
    fillerLength = 5
    alignment = character + "-" * fillerLength + character + "-" * fillerLength + character
    indices = [[0], [fillerLength+1], [len(alignment)-1]]
    differences = [
        (0, 1, character),
        (fillerLength+1, fillerLength+2, character),
        (len(alignment)-1, len(alignment), character),
    ]
    differencesOutput = consensusStrategyPairwiseFunctions.identify_differences_from_indices("mutation", indices, alignment)
    assert differencesOutput == differences

def test_consensus_strategy_pairwise_functions_inject_difference_into_sequence_insertions():
    insert = "-"
    string = "A"*10
    stringWithFrontInsert = insert + string
    frontInsertDifference = (0,0,insert)
    stringWithFrontInsertOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(string, frontInsertDifference)
    assert stringWithFrontInsertOutput == stringWithFrontInsert

    stringWithBackInsert = string + insert
    backInsertDifference = (len(stringWithBackInsert)-1,len(stringWithBackInsert)-1,insert)
    stringWithBackInsertOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(string, backInsertDifference)
    assert stringWithBackInsertOutput == stringWithBackInsert

    halfwayIndex = len(string)//2
    stringWithMiddleInsert = string[:halfwayIndex] + insert + string[halfwayIndex:]
    middleInsertDifference = (halfwayIndex, halfwayIndex, insert)
    stringWithMiddleInsertOutput = consensusStrategyPairwiseFunctions.inject_difference_into_sequence(string, middleInsertDifference)
    assert stringWithMiddleInsertOutput == stringWithMiddleInsert

'''
def test_consensus_strategy_pairwise_functions_apply_difference_into_sequence_deletions():
    insert = "-"
    string = "A"*10
    stringWithFrontInsert = insert + string
    frontDeletionDifference = (0,1,"")
    stringOutput = consensusStrategyPairwiseFunctions.apply_difference_into_sequence(stringWithFrontInsert, frontDeletionDifference)
    assert stringOutput == string

    stringWithBackInsert = string + insert
    backDeletionDifference = (len(stringWithBackInsert)-2,len(stringWithBackInsert)-1,"")
    print(stringWithBackInsert)
    print(string)
    stringOutput = consensusStrategyPairwiseFunctions.apply_difference_into_sequence(stringWithBackInsert, backDeletionDifference)
    assert stringOutput == string

    halfwayIndex = len(string)//2
    stringWithMiddleInsert = string[:halfwayIndex] + insert + string[halfwayIndex:]
    middleDeletionDifference = (halfwayIndex, halfwayIndex+1, "")
    stringOutput = consensusStrategyPairwiseFunctions.apply_difference_into_sequence(stringWithMiddleInsert, middleDeletionDifference)
    assert stringOutput == string
'''
