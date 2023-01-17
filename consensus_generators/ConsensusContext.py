import consensus_generators.ConsensusStrategy as cs

# import ConsensusStrategy as cs


class ConsensusContext:
    def __init__(self, strategy: str):
        self._strategy_types = {
            "pairwise": cs.PairwiseStrategy(),
            "lamassemble": cs.LamassembleStrategy(),
            "medaka": cs.MedakaStrategy(),
        }
        self._strategy = self._strategy_types[strategy]

    @property
    def strategy(self) -> cs.ConsensusStrategy:
        return self._strategy

    @strategy.setter
    def strategy(self, strategy: cs.ConsensusStrategy) -> None:
        self._strategy = strategy

    def generate_consensus_sequences(self, binPaths) -> None:
        return self._strategy.generate_consensus_sequences(binPaths)

    def benchmark_binned_sequences(self, binPath):
        return self._strategy.benchmark_binned_sequences(binPath)
