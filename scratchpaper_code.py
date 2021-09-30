import re


umi_pattern = re.compile('[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}')
example = 'GAGTAAGACAGGTCGGGGCAATGCTTTGAACCAAAGT'
example = 'TGATAGGTCAATTCAGTCCGCCGGATCGAACCAGGGA'
print(not umi_pattern.match(example))
