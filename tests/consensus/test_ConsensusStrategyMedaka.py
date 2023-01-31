from pytestConsensusFixtures import consensusSequence, targetSequences, targetSequenceRecords, simpleInsert, simpleString, middleInsertIndex, targetSequenceDifferences
import pytest
from consensus import ConsensusStrategyMedaka

@pytest.fixture
def consensusStrategyMedaka():
    return ConsensusStrategyMedaka.ConsensusStrategyMedaka()
