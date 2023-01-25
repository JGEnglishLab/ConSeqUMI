from abc import ABC, abstractmethod

class ConsensusStrategy(ABC):

    @abstractmethod
    def generate_consensus_sequence_from_records(self, binRecords: list) -> str:
        pass
