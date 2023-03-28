from abc import ABC, abstractmethod
import random
from Levenshtein import distance
import time
from ConSeqUMI.Printer import Printer
from multiprocessing import Process, Queue

class ConsensusStrategy(ABC):
    @abstractmethod
    def generate_consensus_algorithm_path_header_insert(self) -> str:
        pass

    @abstractmethod
    def generate_consensus_sequence_from_biopython_records(
        self, binRecords: list
    ) -> str:
        pass

    def generate_consensus_algorithm_path_header(self, processName: str):
        return (
            processName
            + "-"
            + self.generate_consensus_algorithm_path_header_insert()
            + time.strftime("-%Y%m%d-%H%M%S")
        )

    def find_consensus_and_add_to_writing_queue(self, queue, randomSampleOfRecords, intervalNumber, iteration,
                                                referenceSequence, numBinRecords):
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
            str(numBinRecords),
        ]
        queue.put(outputList)

    def benchmark_sequence_generator(
        self,
        queue: Queue,
        referenceSequence: str,
        binRecords: list,
        intervalNumbers: list,
        iterations: int,
    ):
        printer = Printer()
        processes = []
        if len(intervalNumbers) == 1:
            intervals = intervalNumbers[0]
            intervalNumbers = [1]
            for i in range(1, len(binRecords) // intervals + 1):
                if i * intervals <= 500:
                    intervalNumbers.append(i * intervals)
        for intervalNumber in intervalNumbers:
            printer(
                f"benchmarking interval: {intervalNumber} ({iterations} iterations)"
            )
            for iteration in range(iterations):
                randomSampleOfRecords = random.sample(binRecords, k=intervalNumber)

                processes.append(Process(target=self.find_consensus_and_add_to_writing_queue, args=(queue, randomSampleOfRecords, intervalNumber, iteration,
                                                referenceSequence, len(binRecords))))
                processes[-1].start()

        for process in processes:
            process.join()