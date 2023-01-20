from general.genomicFunctions import convert_IUPAC_to_regular_expression
from config import UMI_PATTERN
from cutadapt.parser import FrontAdapter, BackAdapter, LinkedAdapter

class UmiExtractor:
    def __init__(self):
        self.umiPattern = f'^{convert_IUPAC_to_regular_expression(UMI_PATTERN)}$'

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
