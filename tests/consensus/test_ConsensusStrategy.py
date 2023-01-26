import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)

from consensus.ConsensusStrategy import ConsensusStrategy

#def test_consensus_strategy_initialization():
