from abc import ABC, abstractmethod
import random
from Levenshtein import distance


class ConsensusStrategy(ABC):
    @abstractmethod
    def generate_consensus_sequence_from_biopython_records(
        self, binRecords: list
    ) -> str:
        pass

    def benchmark_sequence_generator(
        self, referenceSequence: str, binRecords: list, intervals: int, iterations: int
    ):
        intervalNumbers = [1]
        for i in range(1, len(binRecords) // intervals + 1):
            intervalNumbers.append(i * intervals)
        if len(intervalNumbers) > 100:
            intervalNumbers = intervalNumbers[:100]
        for intervalNumber in intervalNumbers:
            for iteration in range(iterations):
                randomSampleOfRecords = random.sample(binRecords, k=intervalNumber)
                if intervalNumber == 1:
                    benchmarkedSequence = str(randomSampleOfRecords[0].seq)
                else:
                    benchmarkedSequence = (
                        self.generate_consensus_sequence_from_biopython_records(
                            randomSampleOfRecords
                        )
                    )
                outputList = [
                    str(intervalNumber),
                    str(iteration),
                    referenceSequence,
                    benchmarkedSequence,
                    str(distance(referenceSequence, benchmarkedSequence)),
                    str(len(binRecords)),
                ]
                yield outputList
