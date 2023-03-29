from ConSeqUMI.consensus.ConsensusStrategyPairwise import (
    ConsensusStrategyPairwise as PairwiseStrategy,
)
from ConSeqUMI.consensus.ConsensusStrategyLamassemble import (
    ConsensusStrategyLamassemble as LamassembleStrategy,
)
from ConSeqUMI.consensus.ConsensusStrategyMedaka import (
    ConsensusStrategyMedaka as MedakaStrategy,
)
from ConSeqUMI.consensus.ConsensusStrategy import ConsensusStrategy as ConsensusStrategy
from concurrent.futures import Future
import typing as T

class ConsensusContext:
    def __init__(self, strategy: str):
        self._strategy_types = {
            "pairwise": PairwiseStrategy(),
            "lamassemble": LamassembleStrategy(),
            "medaka": MedakaStrategy(),
        }
        self._strategy = self._strategy_types[strategy]

    @property
    def strategy(self) -> ConsensusStrategy:
        return self._strategy

    @strategy.setter
    def strategy(self, strategy: ConsensusStrategy) -> None:
        self._strategy = strategy

    def generate_consensus_algorithm_path_header(self, processName: str):
        return self._strategy.generate_consensus_algorithm_path_header(processName)

    def generate_consensus_record_from_biopython_records(
        self, binRecords: list
    ) -> str:
        return self._strategy.generate_consensus_record_from_biopython_records(
            binRecords
        )

    def populate_future_processes_with_benchmark_tasks(
        self, futureProcesses: T.List[Future], numProcesses: int, referenceSequence: str, binRecords: list, intervals: int, iterations: int
    ):
        self._strategy.populate_future_processes_with_benchmark_tasks(
             futureProcesses, numProcesses, referenceSequence, binRecords, intervals, iterations
        )
