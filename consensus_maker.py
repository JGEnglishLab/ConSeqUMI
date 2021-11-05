from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import numpy as np
from collections import Counter
import unittest
from Bio import AlignIO
import subprocess

def remove_chimeras_from_umi_pairs(starcode1Path, starcode2Path, output):
    s1 = pd.read_csv(starcode1Path, sep='\t', header=None)
    s1UMI = s1.iloc[:,0]
    s1Indices = [set([int(y) for y in x.split(',')]) for x in list(s1.iloc[:,2])]
    s2 = pd.read_csv(starcode2Path, sep='\t', header=None)
    s2UMI = s2.iloc[:,0]
    s2Indices = [set([int(y) for y in x.split(',')]) for x in list(s2.iloc[:,2])]

    umi1List = []
    umi2List = []
    indicesList = []
    for i in range(len(s1UMI)):
        umi1 = s1UMI[i]
        indices1 = s1Indices[i]
        for j in range(len(s2UMI)):
            umi2 = s2UMI[j]
            indices2 = s2Indices[j]
            intersect = indices1.intersection(indices2)
            if len(intersect) != 0:
                umi1List.append(umi1)
                umi2List.append(umi2)
                indicesList.append(intersect)

    lengths = [len(i) for i in indicesList]
    lengths, indicesList, umi1List, umi2List = zip(*sorted(zip(lengths, indicesList, umi1List, umi2List), reverse=True))
    umi1Set = set()
    umi2Set = set()
    remove = []
    for i in range(len(indicesList)):
        umi1 = umi1List[i]
        umi2 = umi2List[i]
        if umi1 in umi1Set or umi2 in umi2Set: remove.append(i)
        else: umi1Set.add(umi1); umi2Set.add(umi2)
    indicesList, umi1List, umi2List = [np.delete(np.array(x),(remove)) for x in [indicesList, umi1List, umi2List]]
    data = []
    for i in range(len(indicesList)): data.append([umi1List[i] + umi2List[i], len(indicesList[i]), ','.join([str(x) for x in sorted(indicesList[i])])])
    df = pd.DataFrame(data)
    df.to_csv(output, sep='\t', index=False, header=False)

def find_consensus_seq_from_list(seqs, index, umiBinPath):
    records = (SeqRecord(Seq(seq.strip()), str(index)) for index,seq in enumerate(seqs))
    infile = umiBinPath.split('.')[0] + str(index) + 'unaligned.fa'
    outfile = umiBinPath.split('.')[0] + str(index) + 'aligned.fa'
    SeqIO.write(records, infile, 'fasta')

    process = subprocess.Popen(['/usr/local/anaconda3/envs/longread/bin/muscle',
     '-in', infile,
     '-out', outfile])
    stdout, stderr = process.communicate()

    alignment = AlignIO.read(outfile, format='fasta')
    #for i in range(len(seqs)): alignment.add_sequence(str(i), seqs[i])
    info = AlignInfo.SummaryInfo(alignment)
    return info.dumb_consensus(threshold=0.7)

def find_consensus_sequences_from_umi_bins(umiBinPath, seqPath):
    umiBins = pd.read_csv(umiBinPath, sep='\t')
    umiBins.columns = ['umi','count','indices']
    trimmedBins = umiBins[umiBins['count'] > 10]

    with open(seqPath, 'r') as file: sequences = np.array(file.readlines())

    finalSeqs = []
    for index, row in trimmedBins.iterrows():
        indices = row['indices'].split(',')
        indices = [int(i)-1 for i in indices]
        cluster = sequences[indices]
        consensusSequence = find_consensus_seq_from_list(cluster, index, umiBinPath)
        finalSeqs.append([row['umi'],str(consensusSequence).strip('\n')])

    return pd.DataFrame(finalSeqs, columns=['UMIBin','ConsensusSequence'])
