import pytest
import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from consensus.ReferenceConsensusGenerator import ReferenceConsensusGenerator
import random

def test_reference_consensus_generator_initializer():
    bufferLength = 10
    sequenceWindowLength = 20
    referenceConsensusGenerator = ReferenceConsensusGenerator(bufferLength=bufferLength, sequenceWindowLength=sequenceWindowLength)
    assert referenceConsensusGenerator.bufferLength == bufferLength
    assert referenceConsensusGenerator.sequenceWindowLength == sequenceWindowLength

    defaultBufferLength = 20
    referenceConsensusGenerator = ReferenceConsensusGenerator()
    assert referenceConsensusGenerator.bufferLength == defaultBufferLength

@pytest.fixture
def consensusSequence():
    random.seed(0)
    consensusSequence = "".join(random.choices("ATGC", k=50))
    return consensusSequence

@pytest.fixture
def readSequences(consensusSequence):
    readSequences = []
    insertIndex = 20
    for i in range(1, 6):
        readSequences.append(consensusSequence[:insertIndex] + i*"A" + consensusSequence[insertIndex:])
    deleteIndex = 30
    for i in range(1, 6):
        readSequences.append(consensusSequence[:deleteIndex] + consensusSequence[deleteIndex + i:])
    mutationIndex = 40
    for i in range(5):
        readSequences.append(consensusSequence[:mutationIndex + i] + "R" + consensusSequence[mutationIndex + i + 1:])
    return readSequences

@pytest.fixture
def referenceConsensusGenerator():
    return ReferenceConsensusGenerator()

def custom_justify_sequence_function(sequence, maxLength, bufferLength):
    tempSequence = sequence[:]
    while len(tempSequence) != maxLength: tempSequence += " "
    tempSequence += " " *bufferLength
    return tempSequence

def test_reference_consensus_generator_justify_left_all_string_lengths_with_buffer(readSequences, referenceConsensusGenerator):
    maxLength = 55
    bufferLength = 20
    justifiedSequences = [custom_justify_sequence_function(sequence, maxLength, bufferLength) for sequence in readSequences]
    justifiedSequencesOutput = referenceConsensusGenerator.justify_left_all_string_lengths_with_buffer(readSequences)
    assert justifiedSequencesOutput == justifiedSequences

def test_reference_consensus_generator_initialize_consensus_sequence_front(consensusSequence, readSequences, referenceConsensusGenerator):
    consensusFront = consensusSequence[:10]
    consensusFrontOutput = referenceConsensusGenerator.initialize_consensus_sequence_front(readSequences)
    assert consensusFrontOutput == consensusFront

def test_reference_consensus_generator_find_next_character_in_sequence(consensusSequence, readSequences, referenceConsensusGenerator):
    subSequenceLength = 10
    consensusFront = consensusSequence[:subSequenceLength]
    readSequenceFronts = [sequence[:subSequenceLength+1] for sequence in readSequences]
    nextCharacterInConsensus = consensusSequence[subSequenceLength]
    nextCharacterInConsensusOutput = referenceConsensusGenerator.find_next_character_in_sequence(readSequenceFronts, consensusFront)
    assert nextCharacterInConsensusOutput == nextCharacterInConsensus

def test_reference_consensus_generator_generate_consensus_sequence(consensusSequence, readSequences, referenceConsensusGenerator):
    consensusSequenceOutput = referenceConsensusGenerator.generate_consensus_sequence(readSequences)
    assert consensusSequenceOutput == consensusSequence
