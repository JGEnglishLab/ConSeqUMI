import consensus_generators.ConsensusStrategy
#import ConsensusStrategy

class ConsensusContext:

    def __init__(self, strategy: str):
        self._strategy_types = {
            'pairwise': consensus_generators.ConsensusStrategy.PairwiseStrategy(),
            'lamassemble': consensus_generators.ConsensusStrategy.LamassembleStrategy(),
            #'pairwise': ConsensusStrategy.PairwiseStrategy(),
            #'lamassemble': ConsensusStrategy.LamassembleStrategy(),
            }
        self._strategy = self._strategy_types[strategy]

    @property
    def strategy(self) -> consensus_generators.ConsensusStrategy.ConsensusStrategy:
    #def strategy(self) -> ConsensusStrategy.ConsensusStrategy:
        return self._strategy

    @strategy.setter
    def strategy(self, strategy: consensus_generators.ConsensusStrategy.ConsensusStrategy) -> None:
    #def strategy(self, strategy: ConsensusStrategy.ConsensusStrategy) -> None:
        self._strategy = strategy

    def generate_consensus_sequences(self, binPaths) -> None:
        return self._strategy.generate_consensus_sequences(binPaths)
