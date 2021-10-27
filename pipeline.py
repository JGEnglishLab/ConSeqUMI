#!/usr/bin/python

import sys
import argparse
from os import path, listdir, mkdir
import umi_extractor as ue
import consensus_maker as cm
import subprocess
from timeit import default_timer as timer


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
    files = [args['input']+file for file in listdir(args['input'])]
    print('----> ' + str(round(timer()-startTime, 2)) + ' identifying adapter start indices')
    UMIExtractor.identify_and_set_front_and_reverse_adapter_start_indices_from_file(files)
    print('----> ' + str(round(timer()-startTime, 2)) + ' extracting UMI sequences')
    excludedCount = UMIExtractor.extract_umi_and_sequences_from_files(files, args['output'])
    print('----> ' + str(round(timer()-startTime, 2)) + ' lines excluded: ' + str(excludedCount))

    print('\nStarcode Binning')
    print('---> ' + str(round(timer()-startTime, 2)) + ' begin process')
    process = subprocess.Popen(['../starcode/starcode',
     '-i', args['output']+ 'umi.txt',
     '-o', args['output']+ 'starcode.txt',
     '--seq-id', '-s'])
    stdout, stderr = process.communicate()

    print('\nConsensus Sequence Generation')
    print('----> ' + str(round(timer()-startTime, 2)) + ' obtaining consensus sequences')
    consensusSequences = cm.find_consensus_sequences_from_umi_bins(args['output']+ 'starcode.txt', args['output']+ 'seq.txt')
    print('----> ' + str(round(timer()-startTime, 2)) + ' writing output')
    consensusSequences.to_csv(args['output']+ 'consensus.csv', index=False)


def set_command_line_settings():
    parser = argparse.ArgumentParser(description='Identifying Consensus Sequences from UMI-binnable nanopore reads.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to folder only containing input Nanopore read fastq files.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path for folder output.')
    parser.add_argument('-a', '--adapters', type=str, required=True, help='A text file with f, F, r, R adapters listed. Defaults to: GAGTGTGGCTCTTCGGAT, ATCTCTACGGTGGTCCTAAATAGT, AATGATACGGCGACCACCGAGATC, and CGACATCGAGGTGCCAAAC, respectively.')
    return parser

def check_for_invalid_input(parser, args):

    if not restricted_file(args['input']): parser.error('The -i or --input argument must be an existing directory')
    files = listdir(args['input'])
    for file in files:
        if not restricted_file(file, permittedTypes=['fq', 'fastq']): parser.error('The directory indicated by the -i or --input argument must only contain fastq files (.fq or .fastq)')
    if not restricted_file(args['output']): os.mkdir(args['output'])
    if not restricted_file(args['adapters'], permittedTypes=['txt']): parser.error('The -a or --adapters argument must be an existing file of type text (.txt) format')
    with open(args['adapters'], 'r') as adapterFile: adapters = [adapter.rstrip() for adapter in adapterFile.readlines()]
    if len(adapters) != 4: parser.error('The -a or --adapters argument file must contain 4 adapters.')
    nucleotideCheck = [[characters in ['A','T','G','C'] for characters in adapter] for adapter in adapters]
    if not all(nucleotideCheck): parser.error('The -a or --adapters argument adapters can only contain the nucleotides A,T,G, and C.')

    if args['input'][-1] != '/': args['input'] += '/'
    if args['output'][-1] != '/': args['output'] += '/'


def restricted_file(x, permittedTypes=[]):
    if not len(permittedTypes):
        if not path.isdir(x): return False
        else: return True
    if x.split('.')[-1].lower() not in permittedTypes and not path.isfile(x): return False
    return True


if __name__ == "__main__":
   main()
