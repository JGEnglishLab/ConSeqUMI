from pytestConsensusFixtures import consensusSequence, readSequences, readSequenceRecords, simpleInsert, simpleString, middleInsertIndex, readSequenceDifferences
import pytest
from consensus import ConsensusStrategyMedaka

@pytest.fixture
def consensusStrategyMedaka():
    return ConsensusStrategyMedaka.ConsensusStrategyMedaka()
