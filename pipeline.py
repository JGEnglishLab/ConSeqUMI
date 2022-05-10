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
import random
from Levenshtein import distance
import pandas as pd
from collections import defaultdict, Counter
from os.path import exists
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from statistics import mean, median


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
    print('----> ' + str(round(timer()-startTime, 2)) + ' lines excluded because adapters didn\'t match: ' + str(excludedCount))
    print('\nStarcode Binning')
    print('----> ' + str(round(timer()-startTime, 2)) + ' begin umi1 process')
    run_starcode(args['output']+ 'umi1.txt', args['output']+ 'starcode1.txt')

    print('----> ' + str(round(timer()-startTime, 2)) + ' begin umi2 process')
    run_starcode(args['output']+ 'umi2.txt', args['output']+ 'starcode2.txt')

    print('----> ' + str(round(timer()-startTime, 2)) + ' remove chimeras')
    cm.remove_chimeras_from_umi_pairs(args['output']+ 'starcode1.txt', args['output']+ 'starcode2.txt', args['output']+ 'starcode_without_chimeras.txt')
    print('----> ' + str(round(timer()-startTime, 2)) + ' bin sequences by UMI pair')
    cm.bin_sequences_by_umi_pair(args['output'] + 'seq.fq', args['output']+ 'starcode_without_chimeras.txt')

    if args['benchmarkClusters']:
        print('----> ' + str(round(timer()-startTime, 2)) + ' bootstrapping')
        dfs = []
        binFiles = [args['output']+x for x in os.listdir(args['output']) if re.match('seq_bin\d+\.fq', x)]
        for i in range(len(binFiles)): print(str(i) + ': ' + str(sorted(binFiles)[i]))
        for binFile in sorted(binFiles):
            tempDf = benchmark_binned_sequences(args['output'], binFile, iteration = 10)
            tempDf.to_csv('.'.join(binFile.split('.')[:-1]) + '_benchmark.csv', index = False)
            dfs.append(tempDf)
        df = pd.concat(dfs)
        print('----> ' + str(round(timer()-startTime, 2)) + ' writing benchmarking output')
        df.to_csv(args['oldOutput'] + 'benchmark.csv', index=False)
    else:

        print('\nConsensus Sequence Generation')
        print('----> ' + str(round(timer()-startTime, 2)) + ' obtaining consensus sequences')

        binFiles = sorted([args['output']+x for x in os.listdir(args['output']) if re.match('seq_bin\d+\.fq', x)])
        pattern = '(\d+)'
        records = medaka_pipeline(args['output'], binFiles, pattern)
        print('----> ' + str(round(timer()-startTime, 2)) + ' writing consensus output')
        with open(args['oldOutput'] + 'consensus.fasta', "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
        if args['variants']:
            print('----> ' + str(round(timer()-startTime, 2)) + ' obtaining variant sequences')
            superBinFiles = cluster_consensus_sequences(args['output'], args['oldOutput'] + 'consensus.fasta', binFiles)
            superPattern = '(_super[_\d]+)'
            finalRecords = medaka_pipeline(args['output'], superBinFiles, superPattern)
            with open(args['oldOutput'] + 'variants.fasta', "w") as output_handle:
                SeqIO.write(finalRecords, output_handle, "fasta")
            print('----> ' + str(round(timer()-startTime, 2)) + ' writing variant output')

def medaka_pipeline(outputDir, binFiles, pattern):
    binPattern = "seq_bin" + pattern + "\.fq"
    binPattern = re.compile(binPattern)
    for binFile in binFiles:
        seq = run_medaka_on_file(outputDir, binFile)
        if seq == 'XXXXX': print('empty sequence: ' + binPattern.search(binFile).group(1)); continue
        os.rename(outputDir + 'consensus.fasta', outputDir + 'consensus' + binPattern.search(binFile).group(1) + '.fasta')

    consPattern = re.compile("consensus" + pattern + "\.fasta")

    consensusFiles = [outputDir+x for x in sorted(os.listdir(outputDir)) if re.match('consensus' + pattern + '.fasta',x)]
    records = []
    for file in consensusFiles:
        for record in SeqIO.parse(file, "fasta"): record.id = consPattern.search(file).group(1); records.append(record)
    return records

def adjust_all_string_lengths(strs, buffer_length):
    max_length = len(max(strs, key = len))
    strs = [x.ljust(max_length+buffer_length) for x in strs]
    return strs

def initialize_consensus_output_string(strs, excerpt_length):
    tempList = [x[:excerpt_length] for x in strs]
    final = max(set(tempList), key = tempList.count)
    return final

def find_next_character(subStrings, tempPattern):
    pattern = re.compile(tempPattern + '.')
    baseCountDict = defaultdict(list)
    all_characters = []
    for x in subStrings:
        r = pattern.search(x)
        if r: baseCountDict[r.group(0)[-1]].append(r.start(0)); all_characters.append(r.group(0)[-1])
    c = Counter(all_characters)
    percent_benchmark = .2
    benchmark = int(len(subStrings)*percent_benchmark)
    common = c.most_common()
    if len(common) > 1 and common[1][1] > benchmark:
        min_val = min([median(baseCountDict[x[0]]) for x in common[:2]])
        res = [x[0] for x in common[:2] if median(baseCountDict[x[0]]) == min_val]
        return res[0]
    return c.most_common()[0][0]

def find_consensus(strs, excerpt_length = 10, buffer_length = 20):
    strs = adjust_all_string_lengths(strs, buffer_length)
    final = initialize_consensus_output_string(strs, excerpt_length)
    for i in range(len(strs[0])-buffer_length):
        subStrings = [x[i:i+buffer_length] for x in strs]
        tempPattern = final[-excerpt_length:]
        nex = find_next_character(subStrings, tempPattern)
        final += nex
    return final.strip(' ')

def run_medaka_on_file(outputDir, binFile, bc = False):
    records = [record for record in SeqIO.parse(binFile, "fastq")]
    draftFile = '.'.join(binFile.split('.')[:-1])+'_draft.fq'
    consSequence_count = defaultdict(int)
    returnSequence = 'XXXXX'
    for i in range(len(records)):
        seqStrs = [str(record.seq) for record in records]
        for i in range(len(seqStrs)):
            while len(seqStrs[i]) != 0 and seqStrs[i][-1] == 'A': seqStrs[i] = seqStrs[i][:-1]

        for i in range(len(seqStrs)):
            if len(seqStrs[i]) == 0: return returnSequence
        draftSeq = find_consensus(seqStrs)
        consSeqs = [draftSeq]
        while len(consSeqs) < 5:
            draftSeq = SeqRecord(Seq(consSeqs[-1]),id='draft seq')
            with open(draftFile, "w") as output_handle: SeqIO.write([draftSeq], output_handle, "fasta")


            process = subprocess.Popen(['medaka_consensus',
             '-i', binFile,
             '-d', draftFile,
             '-o', outputDir])
            stdout, stderr = process.communicate()
            print(f'Draft File: {draftFile}')
            print(f'Bin File: {binFile}')
            if exists(outputDir + 'consensus.fasta'): consensusRecords = [record for record in SeqIO.parse(outputDir + 'consensus.fasta', "fasta")]
            else: return returnSequence
            os.remove(outputDir + 'calls_to_draft.bam')
            os.remove(outputDir + 'calls_to_draft.bam.bai')
            os.remove(outputDir + 'consensus_probs.hdf')
            os.remove(draftFile)
            os.remove(draftFile + '.mmi')


            if len(consensusRecords) == 0:
                os.remove(outputDir + 'consensus.fasta')
                return returnSequence

            consensusSequence = str(consensusRecords[0].seq)
            if bc: os.remove(outputDir + 'consensus.fasta')
            if consSeqs[-1] == consensusSequence: return consensusSequence
            consSeqs.append(consensusSequence)
        return consensusSequence
        consSequence_count[consensusSequence] += 1
        if consSequence_count[consensusSequence] >= 3: break

    return consensusSequence

def benchmark_binned_sequences(outDir, binPath, iteration = 100):
    records = [record for record in SeqIO.parse(binPath, "fastq")]
    fullData = []
    if len(records) >= 500:
        tempBinPath = outDir + 'temp_bin.fq'
        referenceSequence = run_medaka_on_file(outDir, binPath, bc = True)

        if referenceSequence == 'XXXXX': return pd.DataFrame(columns = ['binPath','clusterSize','iteration','referenceSequence','tempSequence','levenshteinDistance', 'clusterNum', 'originalClusterSize'])
        data = []
        sequenceData = []
        clusterSizes = [1]
        for i in range(1, len(records)//iteration+1):
            #if i*iteration == 500: clusterSizes.append(i*iteration)
            clusterSizes.append(i*iteration)
        clusterSizes = clusterSizes[:30] #only go to 300
        if len(clusterSizes) > 100: clusterSizes = clusterSizes[100:]
        for i in clusterSizes:
            levenshteinDistances = []
            sequences = []
            for j in range(30):
                tempRecords = random.sample(records, k=i)
                if i == 1: tempSequence = str(tempRecords[0].seq)
                else:
                    with open(tempBinPath, "w") as output_handle:
                        SeqIO.write(tempRecords, output_handle, "fastq")
                    tempSequence = run_medaka_on_file(outDir, tempBinPath, bc = True)
                print('*'*20)
                print('binFile: ' + binPath)
                print('clusterSize: ' + str(i))
                print('iteration: ' + str(j))
                print('distance: ' + str(distance(referenceSequence, tempSequence)))
                print('refSeq: ' + referenceSequence)
                print('temSeq: ' + tempSequence)
                print('*'*20)
                levDist = distance(referenceSequence, tempSequence)
                clusterNum = int(re.search(r'seq_bin(\d+).fq', binPath).group(1))
                originalClusterSize = pd.read_csv(outDir + 'starcode_without_chimeras.txt', sep='\t', header=None).iloc[:,1][clusterNum]

                fullData.append([binPath, i, j, referenceSequence, tempSequence, levDist, clusterNum, originalClusterSize])
                levenshteinDistances.append(distance(referenceSequence, tempSequence))
                sequences.append(tempSequence)

            data.append([i] + levenshteinDistances)
            sequenceData.append([i] + sequences)
    fullDf = pd.DataFrame(fullData, columns = ['binPath','clusterSize','iteration','referenceSequence','tempSequence','levenshteinDistance', 'clusterNum', 'originalClusterSize'])
    return fullDf




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
    parser.add_argument('-bc', '--benchmarkClusters', action="store_true", help='A flag indicating we want to benchmark the optimal cluster size required to generate an accurate consensus sequence.')
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
     '--seq-id',
     #'-s',
     #'-d 5',
     ])
    stdout, stderr = process.communicate()


if __name__ == "__main__":
   main()
