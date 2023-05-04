from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

inputFile = "/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/data/flag_problematic_nucleotides/ConSeqUMI-consensus-20230424-144417/consensus-lamassemble-20230424-144417.fastq"
outputDir = "/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/data/flag_problematic_nucleotides/consensusFiles/"
consRecords = list(SeqIO.parse(inputFile, "fastq"))
countLength = len(str(len(consRecords)))

for i in range(len(consRecords)):
    record = consRecords[i]
    with open(
        outputDir + f"consensus{str(i).zfill(countLength)}.fastq", "w"
    ) as output_handle:
        recordList = [record]
        SeqIO.write(recordList, output_handle, "fastq")
