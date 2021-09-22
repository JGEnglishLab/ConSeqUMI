from Bio.Seq import Seq
from Bio import SeqIO
from cutadapt.parser import AdapterParser
import re

class UMIBinner():

    def __init__(self):
        self.umi_pattern = re.compile('[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}')

    def set_adapters_for_future_matching(self, forwardAdapter1, forwardAdapter2, reverseAdapter1, reverseAdapter2):
        parser = AdapterParser(
                max_errors=0.2, min_overlap=11, read_wildcards=False,
                adapter_wildcards=False, indels=False)

        forward = forwardAdapter1 + '...' + forwardAdapter2
        forward_pair = str(Seq(reverseAdapter2).reverse_complement()) + '...' + str(Seq(reverseAdapter1).reverse_complement())
        reverse = reverseAdapter1 + '...' + reverseAdapter2
        reverse_pair = str(Seq(forwardAdapter2).reverse_complement()) + '...' + str(Seq(forwardAdapter1).reverse_complement())

        self.forward, self.forward_pair, self.reverse, self.reverse_pair = parser.parse_multi(
            [('front',forward), ('back',forward_pair), ('front',reverse), ('back',reverse_pair)]
        )

    def match_sequence_to_adapters(self,seq, strict=False):
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
            if match1 is None or match2 is None: return None

        trimmedRead = match1.trimmed(read1) + match2.trimmed(read2)
        if reverse: trimmedRead = str(Seq(trimmedRead).reverse_complement())
        #if strict and not self.umi_pattern.match(trimmedRead): return None
        return trimmedRead

    def extract_umi_pairs_from_file(self, filePath):
        with open(filePath) as file:
            fileParser = SeqIO.parse(file, 'fastq')
            for record in fileParser:
                seq = str(record.seq)
                seq = self.match_sequence_to_adapters(seq)
                #if seq is not None: yield seq
                yield seq

# NOTE: Porechop and filtlong occur here. It looks like they were not used in the first test. Circle back to later.
