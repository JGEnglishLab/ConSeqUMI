from abc import ABC, abstractmethod
import re
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
from statistics import mean, median
import numpy as np
import subprocess
from io import StringIO
import os
from os.path import exists
import pandas as pd
import random
from Levenshtein import distance


class ConsensusStrategy(ABC):
    def __init__(self):
        self.binPattern = re.compile("seq_bin(\d+)\.fq")

    @abstractmethod
    def generate_consensus_sequence_from_records(self, binRecords: list) -> str:
        pass

    def generate_consensus_sequences(self, binPaths: list):
        self.outputDir = "/".join(binPaths[0].split("/")[:-1]) + "/"
        records = []
        for i in range(len(binPaths)):
            binPath = binPaths[i]
            binNum = self.binPattern.search(binPath).group(1)
            binRecords = [record for record in SeqIO.parse(binPath, "fastq")]
            consensusSeq = self.generate_consensus_sequence_from_records(binRecords)
            seqRecord = SeqRecord(Seq(consensusSeq), id=binNum)
            records.append(seqRecord)
            if i % 1 == 0:
                print(f"{i+1} / {len(binPaths)}", flush=True)
        return records

    def benchmark_binned_sequences(self, binPath, iteration):
        self.outputDir = "/".join(binPath.split("/")[:-1]) + "/"
        binNum = self.binPattern.search(binPath).group(1)
        records = [record for record in SeqIO.parse(binPath, "fastq")]
        fullData = []
        columns = [
            "binPath",
            "clusterSize",
            "iteration",
            "referenceSequence",
            "tempSequence",
            "levenshteinDistance",
            "clusterNum",
            "originalClusterSize",
        ]

        if len(records) >= 300:
            backupFile = self.outputDir + f"backup_benchmark{binNum}.csv"
            with open(backupFile, "w") as f:
                f.write(",".join(columns))
                f.write("\n")

                referenceSequence = self.generate_consensus_sequence_from_records(
                    records
                )
                data = []
                sequenceData = []
                clusterSizes = [1]
                for i in range(1, len(records) // iteration + 1):
                    # if i*iteration == 500: clusterSizes.append(i*iteration)
                    clusterSizes.append(i * iteration)
                if len(clusterSizes) > 100:
                    clusterSizes = clusterSizes[:100]
                for i in clusterSizes:
                    levenshteinDistances = []
                    sequences = []
                    for j in range(100):
                        tempRecords = random.sample(records, k=i)
                        if i == 1:
                            tempSequence = str(tempRecords[0].seq)
                        else:
                            tempSequence = (
                                self.generate_consensus_sequence_from_records(
                                    tempRecords
                                )
                            )
                        print("*" * 20)
                        print("binFile: " + binPath)
                        print("clusterSize: " + str(i))
                        print("iteration: " + str(j))
                        print(
                            "distance: "
                            + str(distance(referenceSequence, tempSequence))
                        )
                        print("refSeq: " + referenceSequence)
                        print("temSeq: " + tempSequence)
                        print("*" * 20, flush=True)
                        levDist = distance(referenceSequence, tempSequence)

                        clusterNum = int(self.binPattern.search(binPath).group(1))
                        originalClusterSize = pd.read_csv(
                            self.outputDir + "starcode_without_chimeras.txt",
                            sep="\t",
                            header=None,
                        ).iloc[:, 1][clusterNum]
                        line = [
                            binPath,
                            i,
                            j,
                            referenceSequence,
                            tempSequence,
                            levDist,
                            clusterNum,
                            originalClusterSize,
                        ]
                        f.write(",".join([str(x) for x in line]))
                        f.write("\n")
                        fullData.append(line)
                        levenshteinDistances.append(
                            distance(referenceSequence, tempSequence)
                        )
                        sequences.append(tempSequence)

                    # if i == 3:
                    #    break
                    data.append([i] + levenshteinDistances)
                    sequenceData.append([i] + sequences)
        fullDf = pd.DataFrame(fullData, columns=columns)
        return fullDf


class ReferenceConsensusGenerator:
    def adjust_all_string_lengths(strs, buffer_length):
        max_length = len(max(strs, key=len))
        strs = [x.ljust(max_length + buffer_length) for x in strs]
        return strs

    def initialize_consensus_output_string(strs, excerpt_length):
        tempList = [x[:excerpt_length] for x in strs]
        final = max(set(tempList), key=tempList.count)
        return final

    def find_next_character(subStrings, tempPattern):
        pattern = re.compile(tempPattern + ".")
        baseCountDict = defaultdict(list)
        all_characters = []
        for x in subStrings:
            r = pattern.search(x)
            if r:
                baseCountDict[r.group(0)[-1]].append(r.start(0))
                all_characters.append(r.group(0)[-1])
        c = Counter(all_characters)
        percent_benchmark = 0.2
        benchmark = int(len(subStrings) * percent_benchmark)
        common = c.most_common()
        if len(common) > 1 and common[1][1] > benchmark:
            min_val = min([median(baseCountDict[x[0]]) for x in common[:2]])
            res = [x[0] for x in common[:2] if median(baseCountDict[x[0]]) == min_val]
            return res[0]
        return c.most_common()[0][0]

    @classmethod
    def find_initial_consensus(self, strs, excerpt_length=10, buffer_length=20):
        strs = self.adjust_all_string_lengths(strs, buffer_length)
        final = self.initialize_consensus_output_string(strs, excerpt_length)
        for i in range(len(strs[0]) - buffer_length):
            subStrings = [x[i : i + buffer_length] for x in strs]
            tempPattern = final[-excerpt_length:]
            nex = self.find_next_character(subStrings, tempPattern)
            final += nex
        return final.strip(" ")


class PairwiseStrategy(ConsensusStrategy):
    def __init__(self):
        super().__init__()

    def find_aligned_differences(self, seq1, seq2):
        aligner = PairwiseAligner()
        alignments = aligner.align(seq1, seq2)
        indices = alignments[0].aligned
        a_idx, b_idx = list(indices[0]), list(indices[1])
        a_idx.insert(0, (0, 0))
        b_idx.insert(0, (0, 0))
        a_idx.append((len(seq1), len(seq1)))
        b_idx.append((len(seq2), len(seq2)))

        diffs = []

        for i in range(len(a_idx) - 1):
            start = a_idx[i][1]
            end = a_idx[i + 1][0]
            insert = seq2[b_idx[i][1] : b_idx[i + 1][0]]
            if start == end and len(insert) == 0:
                continue
            diffs.append((start, end, insert))
        return diffs

    def find_all_diffs_between_candidate_and_binned_seqs(
        self, candidateSeq, binRecords
    ):
        diffs = []
        for binRecord in binRecords:
            tempDiffs = self.find_aligned_differences(candidateSeq, str(binRecord.seq))
            tempDiffs = [(x[0], x[1], x[2], binRecord.id) for x in tempDiffs]
            diffs.extend(tempDiffs)
        return diffs

    def update_candidate_seq_by_common_diffs(self, candidateSeq, binRecords):
        diffs = self.find_all_diffs_between_candidate_and_binned_seqs(
            candidateSeq, binRecords
        )
        finalSeq = candidateSeq[:]

        tempDiffs = [(x[0], x[1], x[2]) for x in diffs]
        most_common_diff = Counter(tempDiffs).most_common(1)[0][0]
        # print(most_common_diff)
        start, end, insert = most_common_diff
        finalSeq = finalSeq[:start] + insert + finalSeq[end:]
        return finalSeq, diffs

    def find_average_aligned_score(self, candidateSeq, binRecords):
        binSeqs = [str(record.seq) for record in binRecords]
        aligner = PairwiseAligner()
        aligner.mismatch_score = -1
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.5
        scores = []
        score_to_ids = []
        for i in range(len(binSeqs)):
            alignments = aligner.align(candidateSeq, binSeqs[i])
            scores.append(alignments[0].score)
            score_to_ids.append([binRecords[i].id, alignments[0].score])
        return mean(scores), score_to_ids

    def generate_consensus_sequence_from_records(self, binRecords):
        diffs = []
        binSeqs = [str(record.seq) for record in binRecords]
        refSeq = ReferenceConsensusGenerator.find_initial_consensus(binSeqs)
        bestScore = -np.inf
        candidateSeq = refSeq[:]
        counter = 0
        curScore, score_to_ids = self.find_average_aligned_score(
            candidateSeq, binRecords
        )
        while curScore >= bestScore:
            bestScore = curScore
            tempSeq, diffs = self.update_candidate_seq_by_common_diffs(
                candidateSeq, binRecords
            )
            curScore, score_to_ids = self.find_average_aligned_score(
                tempSeq, binRecords
            )
            if curScore >= bestScore:
                candidateSeq = tempSeq
            counter += 1
            if counter > 15:
                return candidateSeq
        return candidateSeq


class LamassembleStrategy(ConsensusStrategy):
    def __init__(self):
        super().__init__()

    def generate_consensus_sequence_from_records(self, binRecords):

        tempRecordFile = self.outputDir + "delete.fq"
        with open(tempRecordFile, "w") as output_handle:
            SeqIO.write(binRecords, output_handle, "fastq")

        child = subprocess.Popen(
            [
                "lamassemble",
                "dependencies_download/promethion.mat",
                tempRecordFile,
                "--end",
                "-g60",
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        child_out = child.communicate()[0].decode("utf8")
        seq_ali = list(SeqIO.parse(StringIO(child_out), "fasta"))
        child.stdin.close()
        os.remove(tempRecordFile)

        return str(seq_ali[0].seq).upper()


class MedakaStrategy(ConsensusStrategy):
    def __init__(self):
        super().__init__()

    def generate_consensus_sequence_from_records(self, binRecords):
        binFile = self.outputDir + "delete.fq"
        with open(binFile, "w") as output_handle:
            SeqIO.write(binRecords, output_handle, "fastq")

        draftFile = self.outputDir + "draft.fq"
        consSequence_count = defaultdict(int)
        returnSequence = ""
        for i in range(len(binRecords)):
            seqStrs = [str(record.seq) for record in binRecords]
            for i in range(len(seqStrs)):
                while len(seqStrs[i]) != 0 and seqStrs[i][-1] == "A":
                    seqStrs[i] = seqStrs[i][:-1]

            for i in range(len(seqStrs)):
                if len(seqStrs[i]) == 0:
                    return returnSequence

            draftSeq = ReferenceConsensusGenerator.find_initial_consensus(seqStrs)
            consSeqs = [draftSeq]
            while len(consSeqs) < 5:

                draftSeq = SeqRecord(Seq(consSeqs[-1]), id="draft seq")
                with open(draftFile, "w") as output_handle:
                    SeqIO.write([draftSeq], output_handle, "fasta")

                process = subprocess.Popen(
                    [
                        "medaka_consensus",
                        "-i",
                        binFile,
                        "-d",
                        draftFile,
                        #'-m', 'r941_min_high_g303',
                        "-o",
                        self.outputDir,
                    ]
                )
                stdout, stderr = process.communicate()
                print(f"Draft File: {draftFile}")
                print(f"Bin File: {binFile}")
                if exists(self.outputDir + "consensus.fasta"):
                    consensusRecords = [
                        record
                        for record in SeqIO.parse(
                            self.outputDir + "consensus.fasta", "fasta"
                        )
                    ]
                else:
                    return returnSequence
                if exists(self.outputDir + "calls_to_draft.bam"):
                    os.remove(self.outputDir + "calls_to_draft.bam")
                if exists(self.outputDir + "calls_to_draft.bam.bai"):
                    os.remove(self.outputDir + "calls_to_draft.bam.bai")
                if exists(self.outputDir + "consensus_probs.hdf"):
                    os.remove(self.outputDir + "consensus_probs.hdf")
                if exists(draftFile + ".mmi"):
                    os.remove(draftFile + ".mmi")
                os.remove(self.outputDir + "consensus.fasta")
                os.remove(binFile)
                os.remove(draftFile)

                if len(consensusRecords) == 0:
                    return returnSequence

                consensusSequence = str(consensusRecords[0].seq)
                if consSeqs[-1] == consensusSequence:
                    return consensusSequence
                consSeqs.append(consensusSequence)
            return consensusSequence
