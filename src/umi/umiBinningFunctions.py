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

def identify_reads_that_are_missing_key_values(topUmis, bottomUmis, targetRecords):
    errorMarkers = []
    numErrorTypes = 4
    for i in range(len(topUmis)):
        errors = [0 for _ in range(numErrorTypes)]
        topUmi = topUmis[i]
        bottomUmi = bottomUmis[i]
        targetRecord = targetRecords[i]
        if targetRecord.name=="adapter not found":
            errors[0] = 1
            errorMarkers.append(errors)
            continue
        if topUmi == "": errors[1] = 1
        if bottomUmi == "": errors[2] = 1
        if len(targetRecord) == 0: errors[3] = 1
        errorMarkers.append(errors)
    return errorMarkers

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

def remove_indices_from_related_lists(listOfLists, removeIndices):
    removeIndicesSet = set(removeIndices)
    editedLists = []
    for l in listOfLists:
        editedList = []
        for idx, ele in enumerate(l):
            if idx not in removeIndicesSet:
                editedList.append(ele)
        editedLists.append(editedList)
    return editedLists
    #return [np.delete(np.array(x), (removeIndices)) for x in listOfLists]

def remove_chimeras_from_umi_pairs_and_return_paired_umi_to_read_records_dict(topUmis, bottomUmis, readIndices, chimeraIndices, targetRecords):
    topUmis, bottomUmis, readIndices = remove_indices_from_related_lists([topUmis, bottomUmis, readIndices], chimeraIndices)
    readRecords = []
    for i in range(len(readIndices)):
        binnedRecords = []
        topUmi = topUmis[i]
        bottomUmi = bottomUmis[i]
        binnedIndices = readIndices[i]
        for j in sorted(binnedIndices):
            record = targetRecords[j-1]
            record.description = f"Top UMI: {topUmi}, Bottom UMI: {bottomUmi}; read number: {j}"
            binnedRecords.append(record)
        readRecords.append(binnedRecords)

    pairedUmiToReadRecords = {(topUmis[i],bottomUmis[i]):readRecords[i] for i in range(len(topUmis))}
    return pairedUmiToReadRecords

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
