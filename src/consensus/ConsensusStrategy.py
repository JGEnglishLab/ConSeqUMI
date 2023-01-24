from abc import ABC, abstractmethod

class ConsensusStrategy():
    def __init__(self): pass

    @abstractmethod
    def generate_consensus_sequence_from_records(self, binRecords: list) -> str:
        pass
