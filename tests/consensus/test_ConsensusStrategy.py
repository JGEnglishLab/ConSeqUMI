import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from consensus.ConsensusStrategy import ConsensusStrategy

def test_consensus_strategy_initialization():
    consensusStrategy = ConsensusStrategy()
    assert consensusStrategy.generate_consensus_sequence_from_records

#def test_consensus_strategy_initialization():
