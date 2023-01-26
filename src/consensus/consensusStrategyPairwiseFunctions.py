import more_itertools as mit

iterable = [2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 20]
[list(group) for group in mit.consecutive_groups(iterable)]

def find_in_string_indices_of_character(string, character):
    individialIndices = [index for index, value in enumerate(string) if value == character]
    indicesGroupedConsecutively = [list(group) for group in mit.consecutive_groups(individialIndices)]
    return indicesGroupedConsecutively

def identify_differences_from_indices(type, indices, alignmentSequence=""):
    differences = []
    format_difference_from_index = format_difference_from_index_function_generator(type)
    for index in indices:
        startIndex, endIndex, insert = format_difference_from_index(index, alignmentSequence)
        differences.append((startIndex,endIndex,insert))
    return differences

def format_difference_from_index_function_generator(type):
    match type:
        case "insertion":
            return format_insertion_difference_from_index
        case "deletion":
            return format_deletion_difference_from_index
        case "mutation":
            return format_mutation_difference_from_index
    return "error"

def format_insertion_difference_from_index(index, alignmentSequence):
    startIndex = index[0]
    endIndex = index[0]
    insert = alignmentSequence[index[0]:index[-1]+1]
    return startIndex, endIndex, insert

def format_deletion_difference_from_index(index, emptySequence):
    startIndex = index[0]
    endIndex = index[-1] + 1
    insert = emptySequence
    return startIndex, endIndex, insert

def format_mutation_difference_from_index(index, alignmentSequence):
    startIndex = index[0]
    endIndex = index[-1] + 1
    insert = alignmentSequence[index[0]:index[-1]+1]
    return startIndex, endIndex, insert

def inject_difference_into_sequence(sequence, difference):
    startIndex, endIndex, insert = difference
    alteredSequence = sequence[:startIndex] + insert + sequence[endIndex:]
    return alteredSequence
