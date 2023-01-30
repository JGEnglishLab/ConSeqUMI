from abc import ABC, abstractmethod

class ConsensusStrategy(ABC):

    @abstractmethod
    def generate_consensus_sequence_from_biopython_records(self, binRecords: list) -> str:
        pass

def convert_biopython_fastq_record_to_string(fastqRecord):
    stringOutput = ""
    stringOutput += f"@{fastqRecord.id}" + " <unknown description>\n"
    stringOutput += f"{str(fastqRecord.seq)}" + "\n"
    stringOutput += "+\n"
    qToPhredScores = [str(chr(score+33)) for score in fastqRecord.letter_annotations["phred_quality"]]
    stringOutput += "".join(qToPhredScores) + "\n"
    return stringOutput
