from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import numpy as np
import subprocess
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from Levenshtein import ratio, distance
import matplotlib.pyplot as plt

def remove_chimeras_from_umi_pairs(starcode1Path, starcode2Path, output, tdd = False):

    s1UMI, s1Indices = gather_umis_and_corresponding_indices_from_starcode(starcode1Path, tdd = tdd)
    s2UMI, s2Indices = gather_umis_and_corresponding_indices_from_starcode(starcode2Path, tdd = tdd)

    umiMatch1, umiMatch2, sharedIndices = sort_umi_pairs_by_number_of_matching_indices(s1UMI, s1Indices, s2UMI, s2Indices)
    umiMatch1, umiMatch2, sharedIndices = remove_duplicate_umis_from_pairs(umiMatch1, umiMatch2, sharedIndices)

    data = []
    for i in range(len(sharedIndices)): data.append([umiMatch1[i] + umiMatch2[i], len(sharedIndices[i]), ','.join([str(x) for x in sorted(sharedIndices[i])])])
    df = pd.DataFrame(data)
    df.to_csv(output, sep='\t', index=False, header=False)

def gather_umis_and_corresponding_indices_from_starcode(starcodePath, tdd = False):
    s1 = pd.read_csv(starcodePath, sep='\t', header=None)
    if isinstance(list(s1.iloc[:,2])[0],int): raise Exception('Fewer that 5 UMI clusters found with more than a single sequence')
    s1UMI = s1.iloc[:,0]
    s1Indices = [set([int(y) for y in x.split(',')]) for x in list(s1.iloc[:,2])]
    remove = []
    for i in range(len(s1Indices)):
     if len(s1Indices) < 10: remove.append(i)

    if not tdd:
        s1UMI, s1Indices = [np.delete(np.array(x),(remove)) for x in [s1UMI, s1Indices]]
        if len(s1Indices) < 5: raise Exception('Fewer that 5 UMI clusters found with more than a single sequence')

    return s1UMI, s1Indices

def sort_umi_pairs_by_number_of_matching_indices(s1UMI, s1Indices, s2UMI, s2Indices):
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
    return umi1List, umi2List, indicesList

def remove_duplicate_umis_from_pairs(umi1List, umi2List, indicesList):
    umi1Set = set()
    umi2Set = set()
    remove = []
    for i in range(len(indicesList)):
        umi1 = umi1List[i]
        umi2 = umi2List[i]
        if umi1 in umi1Set or umi2 in umi2Set: remove.append(i)
        else: umi1Set.add(umi1); umi2Set.add(umi2)
    indicesList, umi1List, umi2List = [np.delete(np.array(x),(remove)) for x in [indicesList, umi1List, umi2List]]

    return umi1List, umi2List, indicesList

def bin_sequences_by_umi_pair(seqPath, starcodePath):
    index_recordID = {}
    with open(seqPath) as handle:
        count = 1
        for record in SeqIO.parse(handle, "fastq"): index_recordID[count] = record.id; count += 1

    starcode = pd.read_csv(starcodePath, sep='\t', header=None)
    starcode = starcode[starcode.iloc[:,1] >= 50]
    starcode = list(starcode.iloc[:,2])
    fq = SeqIO.index(seqPath, "fastq")

    for i in range(len(starcode)):
        indices = [int(y) for y in starcode[i].split(',')]
        records = [fq[index_recordID[j]] for j in indices]
        outputPath = '.'.join(seqPath.split('.')[:-1]) + '_bin' + str(i) + '.fq'
        with open(outputPath, "w") as output_handle:
            SeqIO.write(records, output_handle, "fastq")
    fq.close()


def make_hamming_distance_matrix(seqs):
    array = np.array(seqs).reshape(-1,1)
    return pdist(np.array(array), lambda x,y: 1-ratio(x[0],y[0]))

def cluster_longread_consensus_sequences(seqs, threshold = 1/20, dendrogramFile=None):
    dist_matrix = make_hamming_distance_matrix(np.array(seqs))
    link_matrix = linkage(dist_matrix, method = 'centroid')
    labels = fcluster(link_matrix, threshold, criterion='distance')
    if dendrogramFile:
        plt.figure()
        dn = dendrogram(link_matrix)
        plt.savefig(dendrogramFile)
    seqs = np.array(seqs)
    for cluster_id in np.unique(labels):
        yield labels==cluster_id
