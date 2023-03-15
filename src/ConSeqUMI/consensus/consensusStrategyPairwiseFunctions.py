import more_itertools as mit
import functools
import operator


def find_in_string_indices_of_character(string, character):
    individialIndices = [
        index for index, value in enumerate(string) if value == character
    ]
    indicesGroupedConsecutively = [
        list(group) for group in mit.consecutive_groups(individialIndices)
    ]
    return indicesGroupedConsecutively


def inject_difference_into_sequence(sequence, difference):
    startIndex, endIndex, insert = difference
    alteredSequence = sequence[:startIndex] + insert + sequence[endIndex:]
    return alteredSequence


def identify_differences_from_indices(
    type, indices, originalSequenceAlignment, differentSequenceAlignment
):
    differences = []
    format_difference_from_index = format_difference_from_index_function_generator(type)
    insertIndices = find_in_string_indices_of_character(originalSequenceAlignment, "-")
    if insertIndices:
        insertIndices = functools.reduce(operator.iconcat, insertIndices)
    for index in indices:
        numInsertsBeforeIndex = sum(
            indelIndex < index[0] for indelIndex in insertIndices
        )
        startIndex, endIndex, insert = format_difference_from_index(
            index, differentSequenceAlignment
        )
        startIndex, endIndex = (
            startIndex - numInsertsBeforeIndex,
            endIndex - numInsertsBeforeIndex,
        )
        differences.append((startIndex, endIndex, insert))
    return differences


def format_difference_from_index_function_generator(type):
    if type == "insertion":
        return format_insertion_difference_from_index
    if type == "deletion":
        return format_deletion_difference_from_index
    if type == "mutation":
        return format_mutation_difference_from_index
    return "error"


def format_insertion_difference_from_index(index, differentSequenceAlignment):
    startIndex = index[0]
    endIndex = index[0]
    insert = differentSequenceAlignment[index[0] : index[-1] + 1]
    return startIndex, endIndex, insert


def format_deletion_difference_from_index(index, differentSequenceAlignment):
    startIndex = index[0]
    endIndex = index[-1] + 1
    insert = ""
    return startIndex, endIndex, insert


def format_mutation_difference_from_index(index, differentSequenceAlignment):
    startIndex = index[0]
    endIndex = index[-1] + 1
    insert = differentSequenceAlignment[index[0] : index[-1] + 1]
    return startIndex, endIndex, insert
