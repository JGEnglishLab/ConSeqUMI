import sys
sys.path.insert(1, '/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/src')
from consensus.ConsensusStrategyPairwise import ConsensusStrategyPairwise

def test_consensus_strategy_pairwise_initialization():
    pairwise = ConsensusStrategyPairwise()
