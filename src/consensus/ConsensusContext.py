from consensus.ConsensusStrategyPairwise import ConsensusStrategyPairwise as PairwiseStrategy
from consensus.ConsensusStrategyLamassemble import ConsensusStrategyLamassemble as LamassembleStrategy
from consensus.ConsensusStrategy import ConsensusStrategy as ConsensusStrategy

class ConsensusContext:
    def __init__(self, strategy: str):
        self._strategy_types = {
            "pairwise": PairwiseStrategy(),
            "lamassemble": LamassembleStrategy(),
        }
        self._strategy = self._strategy_types[strategy]

    @property
    def strategy(self) -> ConsensusStrategy:
        return self._strategy

    @strategy.setter
    def strategy(self, strategy: ConsensusStrategy) -> None:
        self._strategy = strategy

    def generate_consensus_sequence_from_biopython_records(self, binRecords: list) -> str:
        return self._strategy.generate_consensus_sequence_from_biopython_records(binRecords)