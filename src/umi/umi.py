import os
import subprocess
from io import StringIO
import pandas as pd
from config import SCOMMAND


def main(args):



    binPath = args["output"] + "bins/"
    os.mkdir(binPath)



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

