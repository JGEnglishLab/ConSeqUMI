from pytestConsensusFixtures import consensusSequence, readSequences, readSequenceRecords, simpleInsert, simpleString, middleInsertIndex, readSequenceDifferences
import pytest
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
from consensus import ConsensusStrategyLamassemble

@pytest.fixture
def consensusStrategyLamassemble():
    return ConsensusStrategyLamassemble.ConsensusStrategyLamassemble()

def test__consensus_strategy_lamassemble__generate_consensus_sequence_from_biopython_records(consensusSequence, readSequenceRecords, consensusStrategyLamassemble):
    consensusSequenceOutput = consensusStrategyLamassemble.generate_consensus_sequence_from_biopython_records(readSequenceRecords)
    assert consensusSequenceOutput == consensusSequence
