import pytest
import random
import sys
import os

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
import pytestConsensusFixtures
from consensus.ReferenceConsensusGenerator import ReferenceConsensusGenerator


def test_reference_consensus_generator_initializer():
    bufferLength = 10
    sequenceWindowLength = 20
    referenceConsensusGenerator = ReferenceConsensusGenerator(
        bufferLength=bufferLength, sequenceWindowLength=sequenceWindowLength
    )
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
def referencetargetSequences(referenceConsensusSequence):
    return pytestConsensusFixtures.generate_read_sequences(referenceConsensusSequence)


@pytest.fixture
def referenceConsensusGenerator():
    return ReferenceConsensusGenerator()


def custom_justify_sequence_function(sequence, maxLength, bufferLength):
    tempSequence = sequence[:]
    while len(tempSequence) != maxLength:
        tempSequence += " "
    tempSequence += " " * bufferLength
    return tempSequence


def test_reference_consensus_generator_justify_left_all_string_lengths_with_buffer(
    referencetargetSequences, referenceConsensusGenerator
):
    maxLength = 55
    bufferLength = 20
    justifiedSequences = [
        custom_justify_sequence_function(sequence, maxLength, bufferLength)
        for sequence in referencetargetSequences
    ]
    justifiedSequencesOutput = (
        referenceConsensusGenerator.justify_left_all_string_lengths_with_buffer(
            referencetargetSequences
        )
    )
    assert justifiedSequencesOutput == justifiedSequences


def test_reference_consensus_generator_initialize_consensus_sequence_front(
    referenceConsensusSequence, referencetargetSequences, referenceConsensusGenerator
):
    consensusFront = referenceConsensusSequence[:10]
    consensusFrontOutput = (
        referenceConsensusGenerator.initialize_consensus_sequence_front(
            referencetargetSequences
        )
    )
    assert consensusFrontOutput == consensusFront


def test_reference_consensus_generator_find_next_character_in_sequence(
    referenceConsensusSequence, referencetargetSequences, referenceConsensusGenerator
):
    subSequenceLength = 10
    consensusFront = referenceConsensusSequence[:subSequenceLength]
    targetSequenceFronts = [
        sequence[: subSequenceLength + 1] for sequence in referencetargetSequences
    ]
    nextCharacterInConsensus = referenceConsensusSequence[subSequenceLength]
    nextCharacterInConsensusOutput = (
        referenceConsensusGenerator.find_next_character_in_sequence(
            targetSequenceFronts, consensusFront
        )
    )
    assert nextCharacterInConsensusOutput == nextCharacterInConsensus


def test_reference_consensus_generator_generate_consensus_sequence(
    referenceConsensusSequence, referencetargetSequences, referenceConsensusGenerator
):
    referenceConsensusSequenceOutput = (
        referenceConsensusGenerator.generate_consensus_sequence(
            referencetargetSequences
        )
    )
    assert referenceConsensusSequenceOutput == referenceConsensusSequence
