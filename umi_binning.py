from Bio.Seq import Seq
from Bio import SeqIO
from cutadapt.parser import FrontAdapter, BackAdapter, LinkedAdapter
import re
import numpy as np
from sklearn.cluster import AffinityPropagation
import distance

class UMIBinner():

    def __init__(self):
        self.umi_pattern = re.compile('[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}')
        self.umi_length = 18
        self.is_umi_patterned = True

    def reverse_complement(self, seq): return str(Seq(seq).reverse_complement())

    def make_linked_adapter(self, adapterSeq1, adapterSeq2, name, max_errors=0.2, min_overlap=11, front_required=True, back_required=True):
        adapter1 = FrontAdapter(adapterSeq1, max_errors=max_errors, min_overlap=min_overlap)
        adapter2 = BackAdapter(adapterSeq2, max_errors=max_errors, min_overlap=min_overlap)
        linkedAdapter = LinkedAdapter(adapter1, adapter2, name=name, front_required=front_required, back_required=back_required)
        return linkedAdapter

    def set_adapters_for_future_matching(self, forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2):
        self.forward = self.make_linked_adapter(forwardAdapter1, forwardAdapter2, 'front')
        self.forward_pair = self.make_linked_adapter(self.reverse_complement(reverseAdapter2), self.reverse_complement(reverseAdapter1), 'front_pair')
        self.reverse = self.make_linked_adapter(reverseAdapter1, reverseAdapter2, 'reverse')
        self.reverse_pair = self.make_linked_adapter(self.reverse_complement(forwardAdapter2), self.reverse_complement(forwardAdapter1), 'reverse_pair')

    def extract_umi_from_linked_adapters(self, seq):

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

        trimmedRead = match1.trimmed(read1) + match2.trimmed(read2)
        if reverse: trimmedRead = self.reverse_complement(trimmedRead)
        if self.is_umi_patterned and not self.umi_pattern.match(trimmedRead): return None
        return trimmedRead

    def get_adapter_indices(self, seq):
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

    def set_adapter_indices(self, indices, indel = 10):
        forward_start_index, reverse_start_index = np.median(indices, axis=0)
        self.forward_adapter_start_index = int(forward_start_index - indel)
        self.reverse_adapter_start_index = int(reverse_start_index - indel)
        self.forward_adapter_stop_index = int(forward_start_index + len(self.forward.front_adapter.sequence) + self.umi_length + len(self.forward.back_adapter.sequence) + indel)
        self.reverse_adapter_stop_index = int(reverse_start_index + len(self.reverse.front_adapter.sequence) + self.umi_length + len(self.reverse.back_adapter.sequence) + indel)

        if self.forward_adapter_start_index < 0: self.forward_adapter_start_index = 0
        if self.reverse_adapter_start_index < 0: self.reverse_adapter_start_index = 0
        print(self.forward_adapter_start_index)
        print(self.forward_adapter_stop_index)
        print(self.reverse_adapter_start_index)
        print(self.reverse_adapter_stop_index)


    def find_consensus_sequences(self, umiSeqs):
        # TODO: currently setting as 'exemplar' sequence, in future do consensus sequence
        consensusSeqs = []
        lev_similarity = -1*np.array([[distance.levenshtein(w1,w2) for w1 in umiSeqs] for w2 in umiSeqs])
        affprop = AffinityPropagation(affinity="precomputed", damping=0.5, random_state=None)
        affprop.fit(lev_similarity)
        for cluster_id in np.unique(affprop.labels_):
            exemplar = umiSeqs[affprop.cluster_centers_indices_[cluster_id]]
            consensusSeqs.append(exemplar)
        return consensusSeqs, affprop

    def identify_adapter_start_end_indices(self, files):
        allIndices = []
        for file in files:
            with open(file) as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    indices = self.get_adapter_indices(str(record.seq))
                    if indices != (-1,-1): allIndices.append(indices)
        self.set_adapter_indices(allIndices)

    def identify_umi_sequences(self, files):
        umi_sequences = []
        for file in files:
            with open(file) as handle:
                umi_sequences.extend([self.extract_umi_from_linked_adapters(str(record.seq)) for record in SeqIO.parse(handle, "fastq")])
        umi_sequences = list(filter(None, umi_sequences))
        return umi_sequences


    def identify_consensus_umi_sequences_from_files(self, files):
        self.identify_adapter_start_end_indices(files)
        umi_sequences = self.identify_umi_sequences(files)
        print('seq1: ' + str(len(umi_sequences)))
        umi_sequences = [x for x in umi_sequences if len(x)==36]
        print('seq2: ' + str(len(umi_sequences)))
        import pickle
        pickle.dump(umi_sequences, open('test/data/umi_sequences.p', 'wb'))
        consensus_sequences, consensus = self.find_consensus_sequences(umi_sequences)
        return consensus
