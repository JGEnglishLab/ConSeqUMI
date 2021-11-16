#!/usr/bin/python

import sys
import argparse
import os
import umi_extractor as ue
import consensus_maker as cm
import subprocess
from timeit import default_timer as timer
import re
from Bio import SeqIO
import numpy as np

def main():

    parser = set_command_line_settings()
    args = vars(parser.parse_args())
    check_for_invalid_input(parser, args)
    startTime = timer()
    print('\nUMI Extraction')
    UMIExtractor = ue.UMIExtractor()
    print('----> ' + str(round(timer()-startTime, 2)) + ' setting adapter objects')
    with open(args['adapters'], 'r') as adapterFile: adapters = [adapter.rstrip() for adapter in adapterFile.readlines()]
    UMIExtractor.set_universal_front_and_reverse_linked_adapters(adapters[0], adapters[1], adapters[2], adapters[3])
    files = [args['input']+file for file in os.listdir(args['input'])]
    print('----> ' + str(round(timer()-startTime, 2)) + ' identifying adapter start indices')
    UMIExtractor.identify_and_set_front_and_reverse_adapter_start_indices_from_file(files)
    print('----> ' + str(round(timer()-startTime, 2)) + ' extracting UMI sequences')
    excludedCount = UMIExtractor.extract_umi_and_sequences_from_files(files, args['output'])
    print('----> ' + str(round(timer()-startTime, 2)) + ' lines excluded: ' + str(excludedCount))

    print('\nStarcode Binning')
    print('----> ' + str(round(timer()-startTime, 2)) + ' begin umi1 process')
    run_starcode(args['output']+ 'umi1.txt', args['output']+ 'starcode1.txt')

    print('----> ' + str(round(timer()-startTime, 2)) + ' begin umi2 process')
    run_starcode(args['output']+ 'umi2.txt', args['output']+ 'starcode2.txt')

    print('----> ' + str(round(timer()-startTime, 2)) + ' remove chimeras')
    cm.remove_chimeras_from_umi_pairs(args['output']+ 'starcode1.txt', args['output']+ 'starcode2.txt', args['output']+ 'starcode_without_chimeras.txt')

    print('----> ' + str(round(timer()-startTime, 2)) + ' bin sequences by UMI pair')
    cm.bin_sequences_by_umi_pair(args['output'] + 'seq.fq', args['output']+ 'starcode_without_chimeras.txt')

    print('\nConsensus Sequence Generation')
    print('----> ' + str(round(timer()-startTime, 2)) + ' obtaining consensus sequences')

    binFiles = [args['output']+x for x in os.listdir(args['output']) if re.match('seq_bin\d+\.fq', x)]
    pattern = '(\d+)'
    records = run_medaka(args['output'], binFiles, pattern)
<<<<<<< HEAD
    print('----> ' + str(round(timer()-startTime, 2)) + ' writing consensus output')
    with open(args['oldOutput'] + 'consensus.fasta', "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    if args['variants']:
        print('----> ' + str(round(timer()-startTime, 2)) + ' obtaining variant sequences')
        superBinFiles = cluster_consensus_sequences(args['output'], args['oldOutput'] + 'consensus.fasta', binFiles)
        superPattern = '(_super[_\d]+)'
        finalRecords = run_medaka(args['output'], superBinFiles, superPattern)
        with open(args['oldOutput'] + 'variants.fasta', "w") as output_handle:
            SeqIO.write(finalRecords, output_handle, "fasta")
        print('----> ' + str(round(timer()-startTime, 2)) + ' writing variant output')
=======
    print('----> ' + str(round(timer()-startTime, 2)) + ' writing output')
    '''
    # Code for variant generation. Currently on the backburner, but code may be added back if the project calls for it.
    with open(args['oldOutput'] + 'consensus.fasta', "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
    print('----> ' + str(round(timer()-startTime, 2)) + ' obtaining variant sequences')
    superBinFiles = cluster_consensus_sequences(args['output'], args['oldOutput'] + 'consensus.fasta', binFiles)
    superPattern = '(_super[_\d]+)'
    finalRecords = run_medaka(args['output'], superBinFiles, superPattern)
    with open(args['oldOutput'] + 'variants.fasta', "w") as output_handle:
        SeqIO.write(finalRecords, output_handle, "fasta")
    print('----> ' + str(round(timer()-startTime, 2)) + ' writing output')
    '''
>>>>>>> 2f94b584da980210a1ac2d4bd87c897820b67422

def make_draft_file(binFilePath, draftFilePath):
    top_record = None
    for record in SeqIO.parse(binFilePath, "fastq"):
        if not top_record: top_record = record; continue
        if len(record.seq) > len(top_record.seq): top_record = record
    records = [top_record]
    with open(draftFilePath, "w") as output_handle:
        SeqIO.write(records, output_handle, "fastq")

def run_medaka(outputDir, binFiles, pattern):
    binPattern = "seq_bin" + pattern + "\.fq"
    binPattern = re.compile(binPattern)
    for binFile in binFiles:
        draftFile = binFile.split('.')[0]+'_draft.fq'
        make_draft_file(binFile, draftFile)
        process = subprocess.Popen(['medaka_consensus',
         '-i', binFile,
         '-d', draftFile,
         '-o', outputDir])
        stdout, stderr = process.communicate()
        os.rename(outputDir + 'consensus.fasta', outputDir + 'consensus' + binPattern.search(binFile).group(1) + '.fq')
        os.remove(outputDir + 'calls_to_draft.bam')
        os.remove(outputDir + 'calls_to_draft.bam.bai')
        os.remove(outputDir + 'consensus_probs.hdf')

    consPattern = re.compile("consensus" + pattern + "\.fq")

    consensusFiles = [outputDir+x for x in sorted(os.listdir(outputDir)) if re.match('consensus' + pattern + '.fq',x)]
    records = []
    for file in consensusFiles:
        for record in SeqIO.parse(file, "fasta"): record.id = consPattern.search(file).group(1); records.append(record)
    return records

def cluster_consensus_sequences(outputDir, consensusFile, binFiles):
    binPattern = "seq_bin(\d+)\.fq"
    binPattern = re.compile(binPattern)
    idBinFileDict = {int(binPattern.search(binFile).group(1)):binFile for binFile in binFiles}
    binConsDict = {}
    for record in SeqIO.parse(consensusFile, "fasta"): binConsDict[int(record.id)] = record
    seqs = [str(binConsDict[id].seq) for id in sorted(binConsDict)]
    seqArray = np.array(seqs)
    finalBinFiles = []
    for group in cm.cluster_longread_consensus_sequences(seqs):
        groupIndices = [i for i,x in enumerate(group) if x]
        records = []
        for i in groupIndices:
            records.extend([record for record in SeqIO.parse(idBinFileDict[i], "fastq")])
        superBinFile = outputDir + 'seq_bin_super_' + '_'.join([str(i) for i in groupIndices]) + '.fq'
        with open(superBinFile, "w") as output_handle:
            SeqIO.write(records, output_handle, "fastq")
        finalBinFiles.append(superBinFile)
    return finalBinFiles






def set_command_line_settings():
    parser = argparse.ArgumentParser(description='Identifying Consensus Sequences from UMI-binnable nanopore reads.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to folder only containing input Nanopore read fastq files.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path for folder output. Folder should not currently exist.')
    parser.add_argument('-a', '--adapters', type=str, required=True, help='A text file with f, F, r, R adapters listed. Defaults to: GAGTGTGGCTCTTCGGAT, ATCTCTACGGTGGTCCTAAATAGT, AATGATACGGCGACCACCGAGATC, and CGACATCGAGGTGCCAAAC, respectively.')
    parser.add_argument('-v', '--variants', action="store_true", help='A flag indicating if variants should be deduced from consensus sequences. For example, if consensus sequences 1, 2, and 3 are generated, and sequences 1 and 3 are the same sequence, the variant file will combine them. The variant output would then have 2 sequences.')
    return parser

def check_for_invalid_input(parser, args):

    if not restricted_file_or_directory(args['input']): parser.error('The -i or --input argument must be an existing directory')
    files = os.listdir(args['input'])
    for file in files:
        if not restricted_file_or_directory(file, permittedTypes=['fq', 'fastq']): parser.error('The directory indicated by the -i or --input argument must only contain fastq files (.fq or .fastq)')
    if restricted_file_or_directory(args['output']): parser.error('The -o or --output argument must indicate a directory that does not exist yet')
    if not restricted_file_or_directory(args['adapters'], permittedTypes=['txt']): parser.error('The -a or --adapters argument must be an existing file of type text (.txt) format')
    with open(args['adapters'], 'r') as adapterFile: adapters = [adapter.rstrip() for adapter in adapterFile.readlines()]
    if len(adapters) != 4: parser.error('The -a or --adapters argument file must contain 4 adapters.')
    nucleotideCheck = [[characters in ['A','T','G','C'] for characters in adapter] for adapter in adapters]
    if not all(nucleotideCheck): parser.error('The -a or --adapters argument adapters can only contain the nucleotides A,T,G, and C.')

    if args['input'][-1] != '/': args['input'] += '/'
    if args['output'][-1] != '/': args['output'] += '/'
    newOutput = args['output'] + 'delete/'
    args['oldOutput'] = args['output']
    args['output'] = newOutput
    os.mkdir(args['oldOutput'])
    os.mkdir(args['output'])


def restricted_file_or_directory(x, permittedTypes=[]):
    if not len(permittedTypes):
        if not os.path.isdir(x): return False
        else: return True
    if x.split('.')[-1].lower() not in permittedTypes and not os.path.isfile(x): return False
    return True

def run_starcode(input, output):
    process = subprocess.Popen(['../starcode/starcode',
     '-i', input,
     '-o', output,
     '--seq-id'])
    stdout, stderr = process.communicate()


if __name__ == "__main__":
   main()
