from general import genomicFunctions
from cutadapt.parser import FrontAdapter, BackAdapter, LinkedAdapter

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
        self.umiPattern = f'^{genomicFunctions.convert_IUPAC_to_regular_expression(umiPattern)}$'

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
        self.topLinkedAdapter_reverseComplement = self.create_linked_adapter(
            genomicFunctions.find_reverse_complement(topBackAdapterSeq, isOnlyStandardNucleotide=True),
            genomicFunctions.find_reverse_complement(topFrontAdapterSeq, isOnlyStandardNucleotide=True),
            "top_reverseComplement",
        )
        self.bottomLinkedAdapter_reverseComplement = self.create_linked_adapter(
            genomicFunctions.find_reverse_complement(bottomBackAdapterSeq, isOnlyStandardNucleotide=True),
            genomicFunctions.find_reverse_complement(bottomFrontAdapterSeq, isOnlyStandardNucleotide=True),
            "bottom_reverseComplement",
        )
