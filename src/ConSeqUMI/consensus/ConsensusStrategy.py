from abc import ABC, abstractmethod
import random
from Levenshtein import distance
import time
from ConSeqUMI.Printer import Printer
from concurrent.futures import ProcessPoolExecutor, Future, as_completed
import typing as T
from Bio.SeqRecord import SeqRecord

class ConsensusStrategy(ABC):
    @abstractmethod
    def generate_consensus_algorithm_path_header_insert(self) -> str:
        pass

    @abstractmethod
    def generate_consensus_record_from_biopython_records(
        self, binRecords: list
    ) -> SeqRecord:
        pass

    def generate_consensus_algorithm_path_header(self, processName: str):
        return (
            processName
            + "-"
            + self.generate_consensus_algorithm_path_header_insert()
            + time.strftime("-%Y%m%d-%H%M%S")
        )

    def find_consensus_and_add_to_writing_queue(self, randomSampleOfRecords, intervalNumber, iteration,
                                                referenceSequence, numBinRecords):
        if intervalNumber == 1:
            benchmarkedRecord = randomSampleOfRecords[0]
        else:
            benchmarkedRecord = (
                self.generate_consensus_record_from_biopython_records(
                    randomSampleOfRecords
                )
            )
        benchmarkedSequence = str(benchmarkedRecord.seq)
        outputList = [
            str(intervalNumber),
            str(iteration),
            referenceSequence,
            benchmarkedSequence,
            str(distance(referenceSequence, benchmarkedSequence)),
            str(numBinRecords),
        ]
        return outputList

    def populate_future_processes_with_benchmark_tasks(
        self,
        futureProcesses: T.List[Future],
        processNum: int,
        referenceSequence: str,
        binRecords: list,
        intervalNumbers: list,
        iterations: int,
    ):

        benchmarkGenerationProcessPool: ProcessPoolExecutor = ProcessPoolExecutor(max_workers=processNum)
        if len(intervalNumbers) == 1:
            intervals = intervalNumbers[0]
            intervalNumbers = [1]
            for i in range(1, len(binRecords) // intervals + 1):
                if i * intervals <= 500:
                    intervalNumbers.append(i * intervals)
        for intervalNumber in intervalNumbers:
            for iteration in range(iterations):
                randomSampleOfRecords = random.sample(binRecords, k=intervalNumber)
                futureProcesses.append(
                benchmarkGenerationProcessPool.submit(self.find_consensus_and_add_to_writing_queue, randomSampleOfRecords, intervalNumber, iteration,
                                                referenceSequence, len(binRecords)))
