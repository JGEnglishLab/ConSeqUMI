import os
import subprocess
from io import StringIO
import pandas as pd
from config import SCOMMAND
from umi.UmiExtractor import UmiExtractor
from umi import umiBinningFunctions
from Bio import SeqIO
from Printer import Printer
import pandas as pd

def main(args):
    printer = Printer()
    umiExtractor = UmiExtractor()
    printer("setting top and bottom linked adapters")
    umiExtractor.set_universal_top_and_bottom_linked_adapters(*args["adapters"])
    printer("extract umis and target sequences from all records")
    rawUmisAndTargetSequences = umiExtractor.extract_umis_and_target_sequences_from_all_records(args["input"])
    printer("identify and remove reads that are missing key values")
    errorMarkers = umiBinningFunctions.identify_reads_that_are_missing_key_values(*rawUmisAndTargetSequences)
    errorIndices = [i for i in range(len(errorMarkers)) if 1 in errorMarkers[i]]
    topRawUmis, bottomRawUmis, targetSequences = umiBinningFunctions.remove_indices_from_related_lists(rawUmisAndTargetSequences, errorIndices)
    printer("create 'data_analysis' folder and add dropped read analysis file")
    dataAnalysisPath = args["output"] + "data_analysis/"
    os.mkdir(dataAnalysisPath)
    readErrorDataFrame = pd.DataFrame(errorMarkers, columns = ["Adapter not found", "Top UMI not found", "Bottom UMI not found", "Target Sequence not found"])
    sequenceIds = [sequence.id for sequence in rawUmisAndTargetSequences[2]]
    readErrorDataFrame.insert(0,"Read ID",sequenceIds)
    readErrorDataFrame.to_csv(dataAnalysisPath + "read_error_summary.csv", index=False)
    if len(topRawUmis) == 0: raise RuntimeError("All provided reads were rejected because no UMIs or target sequences were identified. Please see the 'data_analysis/read_error_summary.csv' file in the output for information on why all reads were rejected.")
    printer("run starcode")
    topUmiToReadIndices = starcode(topRawUmis, dataAnalysisPath + "starcode_output_for_top_umis.csv")
    bottomUmiToReadIndices = starcode(bottomRawUmis, dataAnalysisPath + "starcode_output_for_bottom_umis.csv")
    printer("pair top and bottom umi starcode results by matching reads")
    starcodeTopUmis, starcodeBottomUmis, readIndices = umiBinningFunctions.pair_top_and_bottom_umi_by_matching_reads(topUmiToReadIndices, bottomUmiToReadIndices)
    printer("identify and remove chimeras")
    chimeraIndices = umiBinningFunctions.identify_chimera_indices(starcodeTopUmis, starcodeBottomUmis)
    pairedUmiToReadRecords = umiBinningFunctions.remove_chimeras_from_umi_pairs_and_return_paired_umi_to_read_records_dict(starcodeTopUmis, starcodeBottomUmis, readIndices, chimeraIndices, targetSequences)
    printer("add chimera analysis file to 'data_analysis' folder")
    chimeraDataFrame = umiBinningFunctions.compile_chimera_data_analysis_data_frame(starcodeTopUmis, starcodeBottomUmis, readIndices, chimeraIndices)
    chimeraDataFrame.to_csv(dataAnalysisPath + "chimera_summary_of_starcode_matches.csv", index=False)
    printer("create and fill 'bins' folder with target sequences binned by umi pairing")
    binPath = args["output"] + "bins/"
    os.mkdir(binPath)
    count = 0
    countLength = len(str(len(pairedUmiToReadRecords)))
    for umis in sorted(pairedUmiToReadRecords, key=lambda k: len(pairedUmiToReadRecords[k]), reverse=True):
        binnedRecords = pairedUmiToReadRecords[umis]
        with open(binPath + f"targetSequenceBin{str(count).zfill(countLength)}.fastq", "w") as output_handle:
            SeqIO.write(binnedRecords, output_handle, "fastq")
        count += 1
    printer("UMI extraction and binning complete")




def starcode(umis, file = None):
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
    if file: starcodeOutput.to_csv(file, index=False)
    child.stdin.close()
    umiToReadIndicesDict = starcodeOutput.set_index("umi").to_dict()["readIndices"]
    for umi, readIndices in umiToReadIndicesDict.items():
        indices = readIndices.split(",")
        indices = {int(i) for i in indices}
        umiToReadIndicesDict[umi] = indices
    return umiToReadIndicesDict

