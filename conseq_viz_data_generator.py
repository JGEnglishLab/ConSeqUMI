import consensus_generators.ConsensusStrategy as cs
import sys
import re
import os
from Bio import SeqIO
import pandas as pd

outputDir = sys.argv[1]
if outputDir[-1] != '/': outputDir += '/'

consensusFile = [outputDir+x for x in os.listdir(outputDir) if re.match('consensus.+\.fasta', x)][0]
consDict = {int(record.id):str(record.seq).upper() for record in SeqIO.parse(consensusFile, "fasta")}

binFiles = sorted([outputDir + 'delete/' + x for x in os.listdir(outputDir + 'delete/') if re.match('seq_bin\d+\.fq', x)])
binPattern = re.compile('seq_bin(\d+)\.fq')

binDict = {}
for binFile in binFiles:
    binNum = int(binPattern.search(binFile).group(1))
    binRecords = [record for record in SeqIO.parse(binFile, "fastq")]
    binDict[binNum] = binRecords

ps = cs.PairwiseStrategy()
diffs = []
for binNum in sorted(consDict):
    print(f'{binNum} / {len(consDict)}')
    tempDiffs = ps.find_all_diffs_between_candidate_and_binned_seqs(consDict[binNum], binDict[binNum])
    tempDiffs = [(x[0], x[1], x[2], binNum, x[3]) for x in tempDiffs]
    diffs.extend(tempDiffs)

df = pd.DataFrame(diffs, columns=['start','end','insert','binNum','seqID'])
df.to_csv(outputDir + 'data_viz.csv', index=False)
