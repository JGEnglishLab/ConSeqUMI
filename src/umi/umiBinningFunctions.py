import numpy as np
import pandas as pd

def pair_top_and_bottom_umi_by_matching_reads(topUmiToReadIndices, bottomUmiToReadIndices):
    topUmis = []
    bottomUmis = []
    matchingReadIndices = []
    for topUmi, topReadIndices in topUmiToReadIndices.items():
        for bottomUmi, bottomReadIndices in bottomUmiToReadIndices.items():
            intersect = topReadIndices.intersection(bottomReadIndices)
            if len(intersect) == 0: continue
            topUmis.append(topUmi)
            bottomUmis.append(bottomUmi)
            matchingReadIndices.append(intersect)

    lengths = [len(x) for x in matchingReadIndices]
    lengths, topUmis, bottomUmis, matchingReadIndices = zip(*sorted(zip(lengths, topUmis, bottomUmis, matchingReadIndices), reverse=True))
    topUmis, bottomUmis, matchingReadIndices = [list(x) for x in [topUmis, bottomUmis, matchingReadIndices]]
    return topUmis, bottomUmis, matchingReadIndices

def identify_chimera_indices(topUmis, bottomUmis):
    previouslyIdentifiedTopUmis = set()
    previouslyIdentifiedBottomUmis = set()
    chimeraIndices = []
    for i in range(len(topUmis)):
        topUmi = topUmis[i]
        bottomUmi = bottomUmis[i]
        if topUmi in previouslyIdentifiedTopUmis or bottomUmi in previouslyIdentifiedBottomUmis:
            chimeraIndices.append(i)
        else:
            previouslyIdentifiedTopUmis.add(topUmi)
            previouslyIdentifiedBottomUmis.add(bottomUmi)
    return chimeraIndices

def remove_chimeras_from_umi_pairs_and_return_paired_umi_to_read_indices_dict(topUmis, bottomUmis, readIndices, chimeraIndices):
    topUmis, bottomUmis, readIndices = [
        np.delete(np.array(x), (chimeraIndices)) for x in [topUmis, bottomUmis, readIndices]
    ]
    pairedUmiToReadIndices = {topUmis[i]+bottomUmis[i]:readIndices[i] for i in range(len(topUmis))}
    return pairedUmiToReadIndices

def compile_chimera_data_analysis_data_frame(topUmis, bottomUmis, readIndices, chimeraIndices):
    data = []
    for i in range(len(topUmis)):
        isChimera = 0
        if i in chimeraIndices: isChimera = 1
        readIndicesString = "/".join([str(index) for index in sorted(readIndices[i])])
        readIndicesLength = len(readIndices[i])
        data.append([
            topUmis[i],
            bottomUmis[i],
            readIndicesLength,
            readIndicesString,
            isChimera,
        ])
    columns = [
        "top UMI",
        "bottom UMI",
        "Number of Reads",
        "Read Identifiers",
        "Not Chimera",
    ]
    chimeraData = pd.DataFrame(data, columns=columns)
    return chimeraData
