import os
import subprocess
from io import StringIO
import pandas as pd
from config import SCOMMAND
from umi.UmiExtractor import UmiExtractor
from umi import umiBinningFunctions
from Bio import SeqIO

def main(args):

    umiExtractor = UmiExtractor()
    umiExtractor.set_universal_top_and_bottom_linked_adapters(*args["adapters"])
    rawUmisAndTargetSequences = umiExtractor.extract_umis_and_target_sequences_from_all_records(args["input"])
    errorMarkers = umiBinningFunctions.identify_reads_that_are_missing_key_values(*rawUmisAndTargetSequences)
    errorIndices = [i for i in range(len(errorMarkers)) if 1 in errorMarkers[i]]
    topRawUmis, bottomRawUmis, targetSequences = umiBinningFunctions.remove_indices_from_related_lists(rawUmisAndTargetSequences, errorIndices)
    topUmiToReadIndices = starcode(topRawUmis)
    bottomUmiToReadIndices = starcode(bottomRawUmis)
    starcodeTopUmis, starcodeBottomUmis, readIndices = umiBinningFunctions.pair_top_and_bottom_umi_by_matching_reads(topUmiToReadIndices, bottomUmiToReadIndices)
    chimeraIndices = umiBinningFunctions.identify_chimera_indices(starcodeTopUmis, starcodeBottomUmis)
    pairedUmiToReadRecords = umiBinningFunctions.remove_chimeras_from_umi_pairs_and_return_paired_umi_to_read_records_dict(starcodeTopUmis, starcodeBottomUmis, readIndices, chimeraIndices, targetSequences)
    binPath = args["output"] + "bins/"
    os.mkdir(binPath)
    count = 0
    for umis in sorted(pairedUmiToReadRecords, key=lambda k: len(pairedUmiToReadRecords[umis]), reverse=True):
        binnedRecords = pairedUmiToReadRecords[umis]
        with open(binPath + f"targetSequenceBin{count}.fastq", "w") as output_handle:
            SeqIO.write(binnedRecords, output_handle, "fastq")
        count += 1




def starcode(umis):
    umisAsTextFileString = "\n".join(umis)
    child = subprocess.Popen(
        SCOMMAND,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )
    child.stdin.write(umisAsTextFileString.encode())
    child_out = child.communicate()[0].decode("utf8")
    starcodeOutput = pd.read_csv(StringIO(child_out), sep="\t", header=None)
    starcodeOutput.columns = ["umi","count","readIndices"]
    child.stdin.close()
    umiToReadIndicesDict = starcodeOutput.set_index("umi").to_dict()["readIndices"]
    for umi, readIndices in umiToReadIndicesDict.items():
        indices = readIndices.split(",")
        indices = {int(i) for i in indices}
        umiToReadIndicesDict[umi] = indices
    return umiToReadIndicesDict

