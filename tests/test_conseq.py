import pytest
import subprocess
import argparse
import re
from unittest.mock import Mock
from tempfile import TemporaryDirectory, NamedTemporaryFile
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
import conseq
from umi.test_UmiExtractor import exampleForwardRecord, exampleReverseRecord, adapterSequences, topUmi, bottomUmi, targetSequence
from pytestConsensusFixtures import consensusSequence, targetSequences, targetSequenceRecords, simpleInsert

def run_command(commandLine):
    commands = commandLine.split()
    child = subprocess.Popen(commands)
    child.communicate()

@pytest.fixture
def parser():
    def mock_error(message):
        raise argparse.ArgumentTypeError(message)
    argparse.ArgumentParser.error = Mock(side_effect = mock_error)
    parser = conseq.set_command_line_settings()
    return parser

@pytest.fixture
def umiFiles(exampleForwardRecord, exampleReverseRecord, adapterSequences):
    class fileObj:
        def __init__(self):
            self.inputDir = TemporaryDirectory(prefix="conseq_input_test_directory_")
            self.outputDir = TemporaryDirectory(prefix="conseq_output_test_directory_")
            self.adapterFile = NamedTemporaryFile(prefix="conseq_adapter_test_file_", suffix=".txt")
            adapters = [
                adapterSequences["topFrontAdapter"],
                adapterSequences["topBackAdapter"],
                adapterSequences["bottomFrontAdapter"],
                adapterSequences["bottomBackAdapter"],
            ]
            with open(self.adapterFile.name, 'w') as f:
                f.write('\n'.join(adapters))
            self.inputForwardFastq = NamedTemporaryFile(prefix="conseq_adapter_test_forward_file_", dir=self.inputDir.name, suffix=".fastq", delete=False)
            with open(self.inputForwardFastq.name, "w") as output_handle:
                SeqIO.write([exampleForwardRecord], output_handle, "fastq")

            self.inputReverseFastq = NamedTemporaryFile(prefix="conseq_adapter_test_reverse_file_", dir=self.inputDir.name, suffix=".fastq", delete=False)
            with open(self.inputReverseFastq.name, "w") as output_handle:
                SeqIO.write([exampleReverseRecord], output_handle, "fastq")
    return fileObj()

@pytest.fixture
def umiArgs(umiFiles):
    umiArgs = [
        "umi",
        "-i",
        umiFiles.inputDir.name + "/",
        "-o",
        umiFiles.outputDir.name + "/",
        "-a",
        umiFiles.adapterFile.name,
    ]
    return umiArgs

@pytest.fixture
def parsedUmiArgs(parser, umiArgs):
    return vars(parser.parse_args(umiArgs))

@pytest.fixture
def outputDirectoryPattern():
    return r"ConSeqUMI-\d{8}-\d{6}.*"

def test__conseq__set_command_line_settings__succeeds_with_umi_args(parsedUmiArgs): pass

def test__conseq__set_command_line_settings__umi_command_succeeds(parsedUmiArgs, umiArgs, exampleForwardRecord, exampleReverseRecord, outputDirectoryPattern, adapterSequences):
    assert set([parsedUmiArgs["input"][0].seq, parsedUmiArgs["input"][1].seq]) == set([exampleForwardRecord.seq, exampleReverseRecord.seq])
    assert re.match(umiArgs[4] + "/" + outputDirectoryPattern, parsedUmiArgs["output"])
    assert os.path.isdir(parsedUmiArgs["output"])
    adapters = [
        adapterSequences["topFrontAdapter"],
        adapterSequences["topBackAdapter"],
        adapterSequences["bottomFrontAdapter"],
        adapterSequences["bottomBackAdapter"],
    ]
    assert parsedUmiArgs["adapters"] == adapters


def test__conseq__set_command_line_settings__unrecognized_command_fails(parser):
    nonExistentCommand = "nonExistentCommand"
    errorOutput = f"argument command: invalid choice: '{nonExistentCommand}' (choose from 'umi', 'cons', 'benchmark')"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args([nonExistentCommand])

def test__conseq__set_command_line_settings__umi_fails_when_does_not_include_input_directory(parser, umiArgs):
    umiArgsWithoutInput = [umiArgs[0]] + umiArgs[3:]
    errorOutput = "the following arguments are required: -i/--input"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithoutInput)

def test__conseq__set_command_line_settings__umi_fails_when_does_not_include_output_directory(parser, umiArgs):
    umiArgsWithoutInput = umiArgs[:3] + umiArgs[5:]
    errorOutput = "the following arguments are required: -o/--output"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithoutInput)

def test__conseq__set_command_line_settings__umi_fails_when_does_not_include_adapter_file(parser, umiArgs):
    umiArgsWithoutInput = umiArgs[:5]
    errorOutput = "the following arguments are required: -a/--adapter"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithoutInput)

def test__conseq__set_command_line_settings__umi_input_directory_fails_when_does_not_exist(parser, umiArgs):
    falseInputDirectory = "/this/path/does/not/exist"
    umiArgsWithFalseInput = umiArgs[:2] + [falseInputDirectory] + umiArgs[3:]
    errorOutput = "The -i or --input argument must be an existing directory."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithFalseInput)

def test__conseq__set_command_line_settings__umi_input_directory_fails_when_it_is_a_file_not_directory(parser, umiArgs):
    inputFile = NamedTemporaryFile(prefix="conseq_adapter_test_input_file_fail_")
    umiArgsWithInputFile = umiArgs[:2] + [inputFile.name] + umiArgs[3:]
    errorOutput = "The -i or --input argument must be a directory, not a file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithInputFile)

def test__conseq__set_command_line_settings__umi_input_directory_fails_when_empty(parser, umiArgs):
    emptyInputDirectory = TemporaryDirectory(prefix="conseq_input_test_directory_empty_")
    umiArgsWithEmptyInput = umiArgs[:2] + [emptyInputDirectory.name] + umiArgs[3:]
    errorOutput = "The -i or --input argument directory must not be empty."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithEmptyInput)

def test__conseq__set_command_line_settings__umi_input_directory_fails_when_contains_non_fastq_file(parser, umiArgs):
    inputDummyFile = NamedTemporaryFile(prefix="conseq_adapter_test_dummy_", dir=umiArgs[2], suffix=".txt", delete=False)
    errorOutput = f"The -i or --input argument directory must only contain fastq files (.fq or .fastq). Offending file: {inputDummyFile.name.split('/')[-1]}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgs)

def test__conseq__set_command_line_settings__umi_output_directory_parameter_has_program_and_time_tag_when_output_directory_does_not_exist(parser, umiArgs, outputDirectoryPattern):
    outputTag = "outputTag"
    umiArgs[4] = umiArgs[4] + outputTag
    umiArgsWithNewOutputTag = umiArgs
    args = vars(parser.parse_args(umiArgsWithNewOutputTag))
    assert os.path.isdir(args["output"])
    parentDirectoryPath = "/".join(args["output"].split("/")[:-2])
    parentDirectoryContents = os.listdir(parentDirectoryPath)
    assert len(parentDirectoryContents) == 1
    assert re.match(outputTag + "_" + outputDirectoryPattern, parentDirectoryContents[0])

def test__conseq__set_command_line_settings__umi_output_directory_fails_when_parent_directory_does_not_exist(parser, umiArgs):
    falseOutputDirectory = "/this/path/does/not/exist"
    umiArgsWithFalseOutput = umiArgs[:4] + [falseOutputDirectory] + umiArgs[5:]
    errorOutput = "The -o or --output argument directory requires an existing parent directory."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithFalseOutput)

def test__conseq__set_command_line_settings__umi_output_directory_fails_when_it_is_a_file_not_directory(parser, umiArgs):
    outputFile = NamedTemporaryFile(prefix="conseq_adapter_test_output_file_fail_")
    umiArgsWithOutputFile = umiArgs[:4] + [outputFile.name] + umiArgs[5:]
    errorOutput = "The -o or --output argument must be a directory, not a file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithOutputFile)

def test__conseq__set_command_line_settings__umi_adapter_file_fails_when_it_is_not_an_existing_file(parser, umiArgs):
    falseAdapterFile = "/this/path/does/not/exist/adapters.txt"
    umiArgs[6] = falseAdapterFile
    umiArgsWithFalseAdapterFile = umiArgs
    errorOutput = "The -a or --adapters argument must be an existing file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithFalseAdapterFile)

def test__conseq__set_command_line_settings__umi_adapter_file_fails_when_file_is_not_txt_file(parser, umiArgs):
    nonTxtFile = NamedTemporaryFile(prefix="conseq_adapter_test_adapter_file_not_text_fail_", suffix=".nottxt")
    umiArgs[6] = nonTxtFile.name
    umiArgsWithNonTxtAdapterFile = umiArgs
    errorOutput = "The -a or --adapters argument must be a text (.txt) file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithNonTxtAdapterFile)

def test__conseq__set_command_line_settings__umi_adapter_file_fails_when_it_contains_more_than_4_lines(parser, umiArgs, adapterSequences):
    adapterFileWithFiveLines = NamedTemporaryFile(prefix="conseq_adapter_test_adapter_file_extra_line_fail_", suffix=".txt")
    adapters = [
        adapterSequences["topFrontAdapter"],
        adapterSequences["topBackAdapter"],
        adapterSequences["bottomFrontAdapter"],
        adapterSequences["bottomBackAdapter"],
        "AAAAAAAAAAAAA",
    ]
    with open(adapterFileWithFiveLines.name, 'w') as f:
        f.write('\n'.join(adapters))
    umiArgs[6] = adapterFileWithFiveLines.name
    umiArgsWithAdapterFileWithFiveLines = umiArgs
    errorOutput = f"The -a or --adapters argument file must contain exactly 4 adapters. Your file contains: {len(adapters)}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithAdapterFileWithFiveLines)

def test__conseq__set_command_line_settings__umi_adapter_file_fails_when_it_contains_an_adapter_with_non_nucleotide_characters(parser, umiArgs, adapterSequences):
    adapterFileWithNonNucleotideAdapter = NamedTemporaryFile(prefix="conseq_adapter_test_adapter_file_non_nuc_fail_", suffix=".txt")
    adapters = [
        adapterSequences["topFrontAdapter"],
        adapterSequences["topBackAdapter"],
        adapterSequences["bottomFrontAdapter"],
        adapterSequences["bottomBackAdapter"] + "R",
    ]
    with open(adapterFileWithNonNucleotideAdapter.name, 'w') as f:
        f.write('\n'.join(adapters))
    umiArgs[6] = adapterFileWithNonNucleotideAdapter.name
    umiArgsWithAdapterFileWithNonNucleotideAdapter = umiArgs
    errorOutput = "The -a or --adapters argument adapters can only contain the nucleotides A,T,G, and C."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithAdapterFileWithNonNucleotideAdapter)

@pytest.fixture
def consFiles(targetSequenceRecords):
    class fileObj:
        def __init__(self):
            self.parentDir = TemporaryDirectory(prefix="conseq_cons_test_parent_directory_")
            self.inputDir = TemporaryDirectory(prefix="conseq_input_test_directory_", dir=self.parentDir.name)

            self.targetSequenceFastq1 = NamedTemporaryFile(prefix="conseq_test_target_sequence_file1_", dir=self.inputDir.name, suffix=".fastq", delete=False)
            with open(self.targetSequenceFastq1.name, "w") as output_handle:
                SeqIO.write(targetSequenceRecords, output_handle, "fastq")

            self.targetSequenceFastq2 = NamedTemporaryFile(prefix="conseq_test_target_sequence_file2_", dir=self.inputDir.name, suffix=".fastq", delete=False)
            with open(self.targetSequenceFastq2.name, "w") as output_handle:
                SeqIO.write(targetSequenceRecords + targetSequenceRecords, output_handle, "fastq")
    return fileObj()


@pytest.fixture
def consArgs(consFiles):
    consArgs = [
        "cons",
        "-i",
        consFiles.inputDir.name + "/",
    ]
    return consArgs

@pytest.fixture
def parsedConsArgs(parser, consArgs):
    consArgs += ["-min", "10"]
    return vars(parser.parse_args(consArgs))

def test__conseq__set_command_line_settings__cons_succeeds_with_cons_args(parsedConsArgs): pass

def test__conseq__set_command_line_settings__cons_defaults_set_correctly(parser, consArgs):
    args = vars(parser.parse_args(consArgs))
    assert len(args["input"]) == 2
    assert args["consensusAlgorithm"] == "pairwise"
    assert args["minimumReads"] == 50

def test__conseq__set_command_line_settings__cons_succeeds_when_consensusAlgorithm_is_pairwise(parser, consArgs):
    consArgs += ["-c", "pairwise"]
    args = vars(parser.parse_args(consArgs))
    assert args["consensusAlgorithm"] == "pairwise"

def test__conseq__set_command_line_settings__cons_succeeds_when_consensusAlgorithm_is_lamassemble(parser, consArgs):
    consArgs += ["-c", "lamassemble"]
    args = vars(parser.parse_args(consArgs))
    assert args["consensusAlgorithm"] == "lamassemble"

'''
def test__conseq__set_command_line_settings__succeeds_when_consensusAlgorithm_is_medaka(parser, consArgs):
    consArgs += ["-c", "medaka"]
    args = vars(parser.parse_args(consArgs))
    assert args["consensusAlgorithm"] == "medaka"
'''

def test__conseq__set_command_line_settings__cons_fails_when_consensusAlgorithm_is_not_recognized(parser, consArgs):
    errorValue = "unidentified"
    consArgs += ["-c", errorValue]
    errorOutput = f"The -c or --consensusAlgorithm argument must be 'pairwise' or 'lamassemble'. Offending value: {errorValue}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(consArgs)

def test__conseq__set_command_line_settings__cons_accepts_minimumReads(parser, consArgs):
    consArgs += ["-min", "12"]
    args = vars(parser.parse_args(consArgs))
    assert args["minimumReads"] == 12


def test__conseq__set_command_line_settings__cons_fails_when_minimumReads_is_not_an_int(parser, consArgs):
    errorValue = "-10.1"
    consArgs += ["-min", errorValue]
    errorOutput = f"The -min or --minimumReads argument must be an integer. Offending value: {errorValue}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(consArgs)

def test__conseq__set_command_line_settings__cons_fails_when_minimumReads_is_negative(parser, consArgs):
    errorValue = "-10"
    consArgs += ["-min", errorValue]
    errorOutput = f"The -min or --minimumReads argument must be greater than or equal to 1. Offending value: {errorValue}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(consArgs)

@pytest.fixture
def benchmarkFiles(consensusSequence, targetSequenceRecords):
    class fileObj:
        def __init__(self):
            self.parentDir = TemporaryDirectory(prefix="conseq_cons_test_parent_directory_")
            self.inputFile = NamedTemporaryFile(prefix="conseq_test_reference_sequence_file_", dir=self.parentDir.name, suffix=".fastq", delete=False)
            self.outputDir = TemporaryDirectory(prefix="conseq_output_test_directory_", dir=self.parentDir.name)
            self.referenceFile = NamedTemporaryFile(prefix="conseq_test_reference_sequence_file_", dir=self.parentDir.name, suffix=".fasta", delete=False)

            with open(self.inputFile.name, "w") as output_handle:
                SeqIO.write(targetSequenceRecords, output_handle, "fastq")

            with open(self.referenceFile.name, "w") as output_handle:
                consensusSequenceRecord = SeqRecord(Seq(consensusSequence), id="consensus")
                SeqIO.write([consensusSequenceRecord], output_handle, "fasta")
    return fileObj()

@pytest.fixture
def benchmarkArgs(benchmarkFiles):
    benchmarkArgs = [
        "benchmark",
        "-i",
        benchmarkFiles.inputFile.name,
        "-o",
        benchmarkFiles.outputDir.name,
    ]
    return benchmarkArgs

@pytest.fixture
def parsedBenchmarkArgs(parser, benchmarkArgs, benchmarkFiles):
    pass

def test__conseq__set_command_line_settings__benchmark_succeeds_with_benchmark_args(parser, benchmarkArgs):
    parser.parse_args(benchmarkArgs)

def test__conseq__set_command_line_settings__benchmark_defaults_set_correctly(parser, benchmarkArgs):
    args = vars(parser.parse_args(benchmarkArgs))
    assert len(args["input"]) == 14
    assert args["consensusAlgorithm"] == "pairwise"
    assert args["reference"] == ""
    assert args["intervals"] == 10
    assert args["iterations"] == 100

def test__conseq__set_command_line_settings__benchmark_fails_when_does_not_include_input_file(parser, benchmarkArgs):
    benchmarkArgsWithoutInput = [benchmarkArgs[0]] + benchmarkArgs[3:]
    errorOutput = "the following arguments are required: -i/--input"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithoutInput)

def test__conseq__set_command_line_settings__benchmark_input_file_fails_when_it_is_a_directory_not_file(parser, benchmarkArgs):
    inputDirectory = TemporaryDirectory(prefix="conseq_cons_test_parent_directory_")
    benchmarkArgsWithInputDirectory = benchmarkArgs[:2] + [inputDirectory.name] + benchmarkArgs[3:]
    errorOutput = "The -i or --input argument must be a file, not a directory."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithInputDirectory)

def test__conseq__set_command_line_settings__benchmark_input_file_fails_when_does_not_exist(parser, benchmarkArgs):
    falseInputFile = "/this/path/does/not/exist.fastq"
    benchmarkArgsWithFalseInput = benchmarkArgs[:2] + [falseInputFile] + benchmarkArgs[3:]
    errorOutput = "The -i or --input argument must be an existing file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithFalseInput)

def test__conseq__set_command_line_settings__benchmark_input_file_fails_when_contains_non_fastq_file(parser, benchmarkArgs):
    inputTextFile = NamedTemporaryFile(prefix="conseq_adapter_test_dummy_", suffix=".txt")
    benchmarkArgsWithInputTextFile = benchmarkArgs[:2] + [inputTextFile.name] + benchmarkArgs[3:]
    errorOutput = "The -i or --input argument file can only be a fastq file (.fq or .fastq)."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithInputTextFile)

def test__conseq__set_command_line_settings__benchmark_fails_when_does_not_include_output_directory(parser, benchmarkArgs):
    benchmarkArgsWithoutInput = benchmarkArgs[:3] + benchmarkArgs[5:]
    errorOutput = "the following arguments are required: -o/--output"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithoutInput)

def test__conseq__set_command_line_settings__benchmark_output_directory_parameter_has_program_and_time_tag_when_output_directory_does_not_exist(parser, benchmarkArgs, outputDirectoryPattern):
    outputTag = "outputTag"
    benchmarkArgs[4] = benchmarkArgs[4] + "/" + outputTag
    benchmarkArgsWithNewOutputTag = benchmarkArgs
    args = vars(parser.parse_args(benchmarkArgsWithNewOutputTag))
    assert os.path.isdir(args["output"])
    outputWithoutPath = args["output"].split("/")[-2]
    assert re.match(outputTag + "_" + outputDirectoryPattern, outputWithoutPath)

def test__conseq__set_command_line_settings__benchmark_output_directory_fails_when_parent_directory_does_not_exist(parser, benchmarkArgs):
    falseOutputDirectory = "/this/path/does/not/exist"
    benchmarkArgsWithFalseOutput = benchmarkArgs[:4] + [falseOutputDirectory] + benchmarkArgs[5:]
    errorOutput = "The -o or --output argument directory requires an existing parent directory."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithFalseOutput)

def test__conseq__set_command_line_settings__benchmark_output_directory_fails_when_it_is_a_file_not_directory(parser, benchmarkArgs):
    outputFile = NamedTemporaryFile(prefix="conseq_adapter_test_output_file_fail_")
    benchmarkArgsWithOutputFile = benchmarkArgs[:4] + [outputFile.name] + benchmarkArgs[5:]
    errorOutput = "The -o or --output argument must be a directory, not a file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithOutputFile)

def test__conseq__set_command_line_settings__benchmark_succeeds_when_consensusAlgorithm_is_pairwise(parser, benchmarkArgs):
    benchmarkArgs += ["-c", "pairwise"]
    args = vars(parser.parse_args(benchmarkArgs))
    assert args["consensusAlgorithm"] == "pairwise"

def test__conseq__set_command_line_settings__benchmark_succeeds_when_consensusAlgorithm_is_lamassemble(parser, benchmarkArgs):
    benchmarkArgs += ["-c", "lamassemble"]
    args = vars(parser.parse_args(benchmarkArgs))
    assert args["consensusAlgorithm"] == "lamassemble"

def test__conseq__set_command_line_settings__benchmark_fails_when_consensusAlgorithm_is_not_recognized(parser, benchmarkArgs):
    errorValue = "unidentified"
    benchmarkArgs += ["-c", errorValue]
    errorOutput = f"The -c or --consensusAlgorithm argument must be 'pairwise' or 'lamassemble'. Offending value: {errorValue}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgs)

def test__conseq__set_command_line_settings__benchmark_reference_file_succeeds_when_provided(parser, benchmarkArgs, benchmarkFiles):
    benchmarkArgsWithReference = benchmarkArgs + ["-r", benchmarkFiles.referenceFile.name]
    args = parser.parse_args(benchmarkArgsWithReference)

def test__conseq__set_command_line_settings__benchmark_reference_file_fails_when_it_is_a_directory_not_file(parser, benchmarkArgs):
    inputDirectory = TemporaryDirectory(prefix="conseq_cons_test_parent_directory_")
    benchmarkArgsWithInputDirectory = benchmarkArgs + ["-r", inputDirectory.name]
    errorOutput = "The -r or --reference argument must be a file, not a directory."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithInputDirectory)

def test__conseq__set_command_line_settings__benchmark_reference_file_fails_when_does_not_exist(parser, benchmarkArgs):
    falseReferenceFile = "/this/path/does/not/exist.fasta"
    benchmarkArgsWithFalseReference = benchmarkArgs + ["-r", falseReferenceFile]
    errorOutput = "The -r or --reference argument must be an existing file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgsWithFalseReference)

def test__conseq__set_command_line_settings__benchmark_reference_file_fails_when_contains_non_fasta_file(parser, benchmarkArgs):
    referenceFastqFile = NamedTemporaryFile(prefix="conseq_adapter_test_dummy_", suffix=".fastq")
    benchmarkArgsWithReferenceFastqFile = benchmarkArgs + ["-r", referenceFastqFile.name]
    errorOutput = "The -r or --reference argument file can only be a fasta file (.fa or .fasta)."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args( benchmarkArgsWithReferenceFastqFile)

def test__conseq__set_command_line_settings__benchmark_accepts_intervals(parser, benchmarkArgs):
    benchmarkArgs += ["-int", "15"]
    args = vars(parser.parse_args(benchmarkArgs))
    assert args["intervals"] == 15

def test__conseq__set_command_line_settings__cons_fails_when_intervals_is_not_an_int(parser, benchmarkArgs):
    errorValue = "-10.1"
    benchmarkArgs += ["-int", errorValue]
    errorOutput = f"The -int or --intervals argument must be an integer. Offending value: {errorValue}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgs)

def test__conseq__set_command_line_settings__cons_fails_when_intervals_is_negative(parser, benchmarkArgs):
    errorValue = "-10"
    benchmarkArgs += ["-int", errorValue]
    errorOutput = f"The -int or --intervals argument must be greater than or equal to 1. Offending value: {errorValue}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgs)

def test__conseq__set_command_line_settings__benchmark_accepts_iterations(parser, benchmarkArgs):
    benchmarkArgs += ["-iter", "15"]
    args = vars(parser.parse_args(benchmarkArgs))
    assert args["iterations"] == 15

def test__conseq__set_command_line_settings__cons_fails_when_iterations_is_not_an_int(parser, benchmarkArgs):
    errorValue = "-10.1"
    benchmarkArgs += ["-iter", errorValue]
    errorOutput = f"The -iter or --iterations argument must be an integer. Offending value: {errorValue}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgs)

def test__conseq__set_command_line_settings__cons_fails_when_iterations_is_negative(parser, benchmarkArgs):
    errorValue = "-10"
    benchmarkArgs += ["-iter", errorValue]
    errorOutput = f"The -iter or --iterations argument must be greater than or equal to 1. Offending value: {errorValue}"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(benchmarkArgs)


