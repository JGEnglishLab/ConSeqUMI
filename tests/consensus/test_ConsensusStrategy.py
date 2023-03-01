from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)