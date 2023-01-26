import pytest
import sys
import pytestConsensusFixtures
sys.path.insert(1, '../../src')
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
def referenceConsensusSequence():
    random.seed(0)
    referenceConsensusSequence = "".join(random.choices("ATGC", k=50))
    return referenceConsensusSequence

@pytest.fixture
def referenceReadSequences(referenceConsensusSequence):
    return pytestConsensusFixtures.generate_read_sequences(referenceConsensusSequence)

@pytest.fixture
def referenceConsensusGenerator():
    return ReferenceConsensusGenerator()

def custom_justify_sequence_function(sequence, maxLength, bufferLength):
    tempSequence = sequence[:]
    while len(tempSequence) != maxLength: tempSequence += " "
    tempSequence += " " *bufferLength
    return tempSequence

def test_reference_consensus_generator_justify_left_all_string_lengths_with_buffer(referenceReadSequences, referenceConsensusGenerator):
    maxLength = 55
    bufferLength = 20
    justifiedSequences = [custom_justify_sequence_function(sequence, maxLength, bufferLength) for sequence in referenceReadSequences]
    justifiedSequencesOutput = referenceConsensusGenerator.justify_left_all_string_lengths_with_buffer(referenceReadSequences)
    assert justifiedSequencesOutput == justifiedSequences

def test_reference_consensus_generator_initialize_consensus_sequence_front(referenceConsensusSequence, referenceReadSequences, referenceConsensusGenerator):
    consensusFront = referenceConsensusSequence[:10]
    consensusFrontOutput = referenceConsensusGenerator.initialize_consensus_sequence_front(referenceReadSequences)
    assert consensusFrontOutput == consensusFront

def test_reference_consensus_generator_find_next_character_in_sequence(referenceConsensusSequence, referenceReadSequences, referenceConsensusGenerator):
    subSequenceLength = 10
    consensusFront = referenceConsensusSequence[:subSequenceLength]
    readSequenceFronts = [sequence[:subSequenceLength+1] for sequence in referenceReadSequences]
    nextCharacterInConsensus = referenceConsensusSequence[subSequenceLength]
    nextCharacterInConsensusOutput = referenceConsensusGenerator.find_next_character_in_sequence(readSequenceFronts, consensusFront)
    assert nextCharacterInConsensusOutput == nextCharacterInConsensus

def test_reference_consensus_generator_generate_consensus_sequence(referenceConsensusSequence, referenceReadSequences, referenceConsensusGenerator):
    referenceConsensusSequenceOutput = referenceConsensusGenerator.generate_consensus_sequence(referenceReadSequences)
    assert referenceConsensusSequenceOutput == referenceConsensusSequence
