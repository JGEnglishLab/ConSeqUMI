import pytest
from pytestConsensusFixtures import consensusSequence, readSequences, readSequenceRecords, simpleInsert, simpleString, middleInsertIndex, readSequenceDifferences
from consensus import ConsensusStrategyLamassemble

@pytest.fixture
def consensusStrategyLamassemble():
    return ConsensusStrategyLamassemble.ConsensusStrategyLamassemble()

def test__consensus_strategy_lamassemble__generate_consensus_sequence_from_biopython_records(consensusSequence, readSequenceRecords, consensusStrategyLamassemble):
    consensusSequenceOutput = consensusStrategyLamassemble.generate_consensus_sequence_from_biopython_records(readSequenceRecords)
    assert consensusSequenceOutput == consensusSequence
