from general import genomicFunctions
from config import UMI_PATTERN
from cutadapt.parser import FrontAdapter, BackAdapter, LinkedAdapter

class UmiExtractor:
    def __init__(self):
        self.umiPattern = f'^{genomicFunctions.convert_IUPAC_to_regular_expression(UMI_PATTERN)}$'

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
        topAdapterSeq1,
        topAdapterSeq2,
        bottomAdapterSeq1,
        bottomAdapterSeq2,
    ):
        self.topLinkedAdapter = self.create_linked_adapter(
            topAdapterSeq1,
            topAdapterSeq2,
            "top",
        )
        self.bottomLinkedAdapter = self.create_linked_adapter(
            bottomAdapterSeq1,
            bottomAdapterSeq2,
            "bottom",
        )
        self.topLinkedAdapter_reverseComplement = self.create_linked_adapter(
            genomicFunctions.find_reverse_complement(topAdapterSeq2),
            genomicFunctions.find_reverse_complement(topAdapterSeq1),
            "top_reverseComplement",
        )
        self.bottomLinkedAdapter_reverseComplement = self.create_linked_adapter(
            genomicFunctions.find_reverse_complement(bottomAdapterSeq2),
            genomicFunctions.find_reverse_complement(bottomAdapterSeq1),
            "bottom_reverseComplement",
        )
