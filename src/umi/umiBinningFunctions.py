

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
    lengths, matchingReadIndices, topUmis, bottomUmis = zip(*sorted(zip(lengths, matchingReadIndices, topUmis, bottomUmis), reverse=True))
    topUmis, bottomUmis, matchingReadIndices = [list(x) for x in [topUmis, bottomUmis, matchingReadIndices]]
    return topUmis, bottomUmis, matchingReadIndices
