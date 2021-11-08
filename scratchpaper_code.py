from Bio import SeqIO

records = []
handle = 'test/data/test_reads.fq'
for record in SeqIO.parse(handle, "fastq"):
    #print(record.id)
    if record.id == '66ca31e2-e5b3-49fd-9d86-8b1824c6d294':
        #print(record.features)
        records.append(record[1:])
        records.append(record.reverse_complement())

with open("test/data/test_reads_subset.fq", "w") as output_handle:
    SeqIO.write(records, output_handle, "fastq")
