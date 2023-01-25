import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from consensus import consensusStrategyPairwiseFunctions


def test_consensus_strategy_pairwise_functions_find_in_string_indices_of_character_for_single_instance_of_character():
    insert = "-"
    string = "A"*10
    stringWithFrontInsert = insert + string
    frontInsertIndices = [[0]]
    frontInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithFrontInsert, character=insert)
    assert frontInsertIndicesOutput == frontInsertIndices

    stringWithBackInsert = string + insert
    backInsertIndices = [[len(stringWithBackInsert)-1]]
    backInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithBackInsert, character=insert)
    assert backInsertIndicesOutput == backInsertIndices

    halfwayIndex = len(string)//2
    stringWithMiddleInsert = string[:halfwayIndex] + insert + string[halfwayIndex:]
    middleInsertIndices = [[len(string)//2]]
    middleInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithMiddleInsert, character=insert)
    assert middleInsertIndicesOutput == middleInsertIndices

def test_consensus_strategy_pairwise_functions_find_in_string_indices_of_character_for_multiple_instances_of_character():
    insert = "-"
    string = "A"*10
    stringWithFrontInsert = insert * 2 + string
    frontInsertIndices = [[0,1]]
    frontInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithFrontInsert, character=insert)
    assert frontInsertIndicesOutput == frontInsertIndices

    stringWithBackInsert = string + insert * 2
    backInsertIndices = [[len(stringWithBackInsert)-2, len(stringWithBackInsert)-1]]
    backInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithBackInsert, character=insert)
    assert backInsertIndicesOutput == backInsertIndices

    halfwayIndex = len(string)//2
    stringWithMiddleInsert = string[:halfwayIndex] + insert * 2 + string[halfwayIndex:]
    middleInsertIndices = [[halfwayIndex, halfwayIndex+1]]
    middleInsertIndicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithMiddleInsert, character=insert)
    assert middleInsertIndicesOutput == middleInsertIndices

def test_consensus_strategy_pairwise_functions_find_in_string_indices_of_character_with_separate_single_and_multiple_instances_of_character():
    insert = "-"
    string = "A"*10
    stringWithDoubleFrontInsert = insert * 2 + string
    stringWithDoubleFrontInsertAndSingleBackInsert = stringWithDoubleFrontInsert + insert
    indices = [[0,1], [len(stringWithDoubleFrontInsertAndSingleBackInsert)-1]]
    indicesOutput = consensusStrategyPairwiseFunctions.find_in_string_indices_of_character(string=stringWithDoubleFrontInsertAndSingleBackInsert, character=insert)
    assert indicesOutput == indices

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
