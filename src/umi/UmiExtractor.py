from umi import umiExtractionFunctions
from cutadapt.parser import FrontAdapter, BackAdapter, LinkedAdapter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class UmiExtractor:
    def __init__(self, *args, **kwargs):
        self.set_umi_pattern(kwargs.get('umiPattern',''))
        if {'topFrontAdapter','topBackAdapter','bottomFrontAdapter','bottomBackAdapter'}.issubset(set(kwargs)):
            self.set_universal_top_and_bottom_linked_adapters(
                kwargs.get('topFrontAdapter',''),
                kwargs.get('topBackAdapter',''),
                kwargs.get('bottomFrontAdapter',''),
                kwargs.get('bottomBackAdapter',''),
            )

    def set_umi_pattern(self, umiPattern):
        self.umiPattern = f'^{umiExtractionFunctions.convert_IUPAC_to_regular_expression(umiPattern)}$'

    def create_linked_adapter(self, sequence1, sequence2, name):
        adapter1 = FrontAdapter(sequence1)
        adapter2 = BackAdapter(sequence2)
        linkedAdapter = LinkedAdapter(
            adapter1,
            adapter2,
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

        self.topLinkedAdapter = self.create_linked_adapter(
            topFrontAdapterSeq,
            topBackAdapterSeq,
            "top",
        )
        self.bottomLinkedAdapter = self.create_linked_adapter(
            bottomFrontAdapterSeq,
            bottomBackAdapterSeq,
            "bottom",
        )

    def find_matches_of_adapters_in_sequence(self, sequence):
        topSequence, bottomSequence = umiExtractionFunctions.extract_top_and_bottom_of_sequence(sequence)
        topMatch = self.topLinkedAdapter.match_to(topSequence)
        bottomMatch = self.bottomLinkedAdapter.match_to(bottomSequence)
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
            return "","",SeqRecord(Seq(""),id="adapter not found")

        topUmi = topMatch.trimmed(topSequence)
        bottomUmi = bottomMatch.trimmed(bottomSequence)
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
