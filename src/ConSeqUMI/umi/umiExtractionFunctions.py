from Bio.Seq import Seq
from cutadapt.adapters import LinkedMatch


def convert_IUPAC_to_regular_expression(nucleotides):
    iupacToRegExDict = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "W": "[AT]",
        "S": "[CG]",
        "R": "[AG]",
        "Y": "[CT]",
        "K": "[GT]",
        "M": "[AC]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
        "N": "[ACGT]",
    }
    regEx = ""
    for n in nucleotides:
        regEx += iupacToRegExDict[n]
    return regEx

def is_IUPAC_nucleotide(sequence):
    allIUPAC_nucleotides = "ACGTWSRYKMBDHVN"
    return set(sequence).issubset([*allIUPAC_nucleotides])


def find_reverse_complement(sequence):
    sequence = sequence.upper()
    if not is_IUPAC_nucleotide(sequence):
        raise ValueError("Provided sequence contains non-nucleotide characters")
    return str(Seq(sequence).reverse_complement())


def extract_top_and_bottom_of_sequence(sequence):
    seqExtractionLength = 200
    topSequence = sequence[:seqExtractionLength]
    bottomSequence_reverseComplement = find_reverse_complement(
        sequence[-seqExtractionLength:]
    )
    return topSequence, bottomSequence_reverseComplement


def find_index_at_end_of_back_adapter(match):
    if match.front_match:
        return match.front_match.rstop + match.back_match.rstop
    else:
        return match.back_match.rstop
