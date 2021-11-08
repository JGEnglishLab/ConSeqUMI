from Bio.Seq import Seq
from cutadapt.parser import FrontAdapter, BackAdapter, LinkedAdapter
import re
import numpy as np
from Bio import SeqIO




class UMIExtractor():

    def __init__(self):
        self.umi_pattern = re.compile('^[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}$')
        self.umi_length = 18
        self.is_umi_patterned = True

    def reverse_complement(self, seq): return str(Seq(seq).reverse_complement())

    def make_linked_adapter(self, adapterSeq1, adapterSeq2, name, max_errors=0.2, min_overlap=11, front_required=True, back_required=True):
        adapter1 = FrontAdapter(adapterSeq1, max_errors=max_errors, min_overlap=min_overlap)
        adapter2 = BackAdapter(adapterSeq2, max_errors=max_errors, min_overlap=min_overlap)
        linkedAdapter = LinkedAdapter(adapter1, adapter2, name=name, front_required=front_required, back_required=back_required)
        return linkedAdapter

    def set_universal_front_and_reverse_linked_adapters(self, forwardAdapterSeq1, forwardAdapterSeq2, reverseAdapterSeq1, reverseAdapterSeq2):
        self.forward = self.make_linked_adapter(forwardAdapterSeq1, forwardAdapterSeq2, 'front')
        self.forward_pair = self.make_linked_adapter(self.reverse_complement(reverseAdapterSeq2), self.reverse_complement(reverseAdapterSeq1), 'front_pair')
        self.reverse = self.make_linked_adapter(reverseAdapterSeq1, reverseAdapterSeq2, 'reverse')
        self.reverse_pair = self.make_linked_adapter(self.reverse_complement(forwardAdapterSeq2), self.reverse_complement(forwardAdapterSeq1), 'reverse_pair')

    def identify_and_set_front_and_reverse_adapter_start_indices_from_file(self, files):
        allIndices = []
        for file in files:
            with open(file) as handle:

                for record in SeqIO.parse(handle, "fastq"):
                    indices = self.get_front_and_reverse_adapter_start_indices(str(record.seq))
                    if indices != (-1,-1): allIndices.append(indices)
        self.set_universal_forward_reverse_adapter_indices(allIndices)

    def get_front_and_reverse_adapter_start_indices(self, seq):
        baseEndNumber = 200
        read1 = seq[:baseEndNumber]
        read2 = seq[-baseEndNumber:]
        reverse = False

        match1 = self.forward.match_to(read1)
        match2 = self.forward_pair.match_to(read2)

        if match1 is None or match2 is None:
            reverse = True
            match1 = self.reverse.match_to(read1)
            match2 = self.reverse_pair.match_to(read2)

        if match1 is None or match2 is None: return -1, -1
        if not reverse: return match1.front_match.rstart, baseEndNumber - match2.front_match.rstop - match2.back_match.rstop
        else: return baseEndNumber - match2.front_match.rstop - match2.back_match.rstop, match1.front_match.rstart

    def set_universal_forward_reverse_adapter_indices(self, indices, indel = 10):
        if len(indices) < 5: raise Exception('Adapter matched to fewer than 5 given sequences')
        forward_start_index, reverse_start_index = np.median(indices, axis=0)
        self.forward_adapter_start_index = int(forward_start_index - indel)
        self.reverse_adapter_start_index = int(reverse_start_index - indel)
        self.forward_adapter_stop_index = int(forward_start_index + len(self.forward.front_adapter.sequence) + self.umi_length + len(self.forward.back_adapter.sequence) + indel)
        self.reverse_adapter_stop_index = int(reverse_start_index + len(self.reverse.front_adapter.sequence) + self.umi_length + len(self.reverse.back_adapter.sequence) + indel)

        if self.forward_adapter_start_index < 0: self.forward_adapter_start_index = 0
        if self.reverse_adapter_start_index < 0: self.reverse_adapter_start_index = 0

    def extract_umi_and_sequence_from_linked_adapters(self, record):
        seq = str(record.seq)
        read1 = seq[self.forward_adapter_start_index:self.forward_adapter_stop_index]
        read2 = seq[-self.reverse_adapter_stop_index:-self.reverse_adapter_start_index]
        if self.reverse_adapter_start_index==0: read2 = seq[-self.reverse_adapter_stop_index:]
        reverse = False

        match1 = self.forward.match_to(read1)
        match2 = self.forward_pair.match_to(read2)

        if match1 is None or match2 is None:
            reverse = True

            read1 = seq[self.reverse_adapter_start_index:self.reverse_adapter_stop_index]
            read2 = seq[-self.forward_adapter_stop_index:-self.forward_adapter_start_index]
            if self.forward_adapter_start_index==0: read2 = seq[-self.forward_adapter_stop_index:]

            match1 = self.reverse.match_to(read1)
            match2 = self.reverse_pair.match_to(read2)
            if match1 is None or match2 is None: return None

        umi1 = match1.trimmed(read1)
        umi2 = match2.trimmed(read2)
        if (len(umi1) > self.umi_length + 2) or (len(umi1) < self.umi_length - 2) or (len(umi2) > self.umi_length + 2) or (len(umi2) < self.umi_length - 2): return None

        if not reverse:
            startIndex = self.forward_adapter_start_index + match1.front_match.rstop + match1.back_match.rstop
            endIndex = self.reverse_adapter_stop_index - match2.front_match.rstart
            centerSeq = record[startIndex:-endIndex]

        else:
            umi1_temp = self.reverse_complement(umi2)
            umi2 = self.reverse_complement(umi1)
            umi1 = umi1_temp
            startIndex = self.reverse_adapter_start_index + match1.front_match.rstop + match1.back_match.rstop
            endIndex = self.forward_adapter_stop_index - match2.front_match.rstart
            centerSeq = record[startIndex:-endIndex]
            centerSeq = centerSeq.reverse_complement(id=True, name=True, description=True)

        return umi1, umi2, centerSeq

    def extract_umi_and_sequences_from_files(self, files, outputDir):
        sequences = []
        for file in files:
            with open(file) as handle:
                sequences.extend([self.extract_umi_and_sequence_from_linked_adapters(record) for record in SeqIO.parse(handle, "fastq")])
        count = len(sequences)
        sequences = list(filter(None, sequences))
        centers = []
        with open(outputDir + "/umi1.txt", "w") as umiFile1, open(outputDir + "/umi2.txt", "w") as umiFile2:
            for umi_center in sequences:
                umi1 = umi_center[0]
                umi2 = umi_center[1]
                center = umi_center[2]
                umiFile1.write(umi1 + '\n')
                umiFile2.write(umi2 + '\n')

                centers.append(center)

        with open(outputDir + "/seq.fq", "w") as centerFile:
                SeqIO.write(centers, centerFile, "fastq")

        return count - len(sequences)
