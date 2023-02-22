from umi import umiExtractionFunctions
from cutadapt.parser import FrontAdapter, BackAdapter, LinkedAdapter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class UmiExtractor:
    def __init__(self, *args, **kwargs):
        if {'topFrontAdapter','topBackAdapter','bottomFrontAdapter','bottomBackAdapter'}.issubset(set(kwargs)):
            self.set_universal_top_and_bottom_linked_adapters(
                kwargs.get('topFrontAdapter',''),
                kwargs.get('topBackAdapter',''),
                kwargs.get('bottomFrontAdapter',''),
                kwargs.get('bottomBackAdapter',''),
            )

    def create_linked_adapter(self, frontAdapterSequence, backAdapterSequence, name):
        frontAdapter = FrontAdapter(frontAdapterSequence, max_errors=0.2, min_overlap=11)
        backAdapter = BackAdapter(backAdapterSequence, max_errors=0.2, min_overlap=11)
        linkedAdapter = LinkedAdapter(
            frontAdapter,
            backAdapter,
            name=name,
            front_required=True,
            back_required=True,
        )
        return linkedAdapter

    def set_universal_top_and_bottom_linked_adapters(
        self,
        topFrontAdapterSeq,
        topBackAdapterSeq,
        bottomFrontAdapterSeq,
        bottomBackAdapterSeq,
    ):
        for seq in [topFrontAdapterSeq,topBackAdapterSeq,bottomFrontAdapterSeq,bottomBackAdapterSeq]:
            if not umiExtractionFunctions.is_only_standard_nucleotide(seq):
                raise ValueError("Provided Adapter Sequences must only contain standard nucleotides (A, T, C, G)")

        self.topAdapter = self.create_linked_adapter(
            topFrontAdapterSeq,
            topBackAdapterSeq,
            "top",
        )
        self.bottomAdapter = self.create_linked_adapter(
            bottomFrontAdapterSeq,
            bottomBackAdapterSeq,
            "bottom",
        )

    def find_matches_of_adapters_in_sequence(self, sequence):
        topSequence, bottomSequence = umiExtractionFunctions.extract_top_and_bottom_of_sequence(sequence)
        topMatch = self.topAdapter.match_to(topSequence)
        bottomMatch = self.bottomAdapter.match_to(bottomSequence)
        return topMatch, bottomMatch

    def extract_umis_and_target_sequence_from_read(self, record):
        tempRecord = record[:]
        sequence = str(record.seq)
        topSequence, bottomSequence = umiExtractionFunctions.extract_top_and_bottom_of_sequence(sequence)
        topMatch, bottomMatch = self.find_matches_of_adapters_in_sequence(sequence)

        if topMatch is None or bottomMatch is None:
            tempRecord = record.reverse_complement()
            tempRecord.id = record.id
            sequence_reverseComplement = umiExtractionFunctions.find_reverse_complement(sequence)
            topSequence, bottomSequence = umiExtractionFunctions.extract_top_and_bottom_of_sequence(sequence_reverseComplement)
            topMatch, bottomMatch = self.find_matches_of_adapters_in_sequence(sequence_reverseComplement)

        if topMatch is None or bottomMatch is None:
            return "","",SeqRecord(Seq(""),name="adapter not found", id=record.id)

        topUmi = umiExtractionFunctions.extract_previously_identified_umi_from_read(topMatch, topSequence)
        bottomUmi = umiExtractionFunctions.extract_previously_identified_umi_from_read(bottomMatch, bottomSequence)
        targetSeqStartIndex = topMatch.front_match.rstop + topMatch.back_match.rstop
        targetSeqEndIndex = bottomMatch.front_match.rstop + bottomMatch.back_match.rstop
        targetSequenceRecord = tempRecord[targetSeqStartIndex:-targetSeqEndIndex]
        return topUmi, bottomUmi, targetSequenceRecord

    def extract_umis_and_target_sequences_from_all_records(self, records):
        topUmis, bottomUmis, targetSequenceRecords = [], [], []
        for record in records:
            topUmi, bottomUmi, targetSequenceRecord = self.extract_umis_and_target_sequence_from_read(record)
            topUmis.append(topUmi)
            bottomUmis.append(bottomUmi)
            targetSequenceRecords.append(targetSequenceRecord)
        return topUmis, bottomUmis, targetSequenceRecords
