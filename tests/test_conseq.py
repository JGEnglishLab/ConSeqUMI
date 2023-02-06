import pytest
import subprocess
import argparse
import re
from unittest.mock import Mock
from tempfile import TemporaryDirectory, NamedTemporaryFile
from Bio import SeqIO
import os

import test_
import conseq
from umi.test_UmiExtractor import exampleForwardRecord, exampleReverseRecord, adapterSequences, topUmi, bottomUmi, targetSequence

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
def files(exampleForwardRecord, exampleReverseRecord, adapterSequences):
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
def umiArgs(files):
    umiArgs = [
        "umi",
        "-i",
        files.inputDir.name + "/",
        "-o",
        files.outputDir.name + "/",
        "-a",
        files.adapterFile.name,
    ]
    return umiArgs

@pytest.fixture
def parsedUmiArgs(parser, umiArgs):
    return vars(parser.parse_args(umiArgs))

@pytest.fixture
def outputDirectoryPattern():
    return r"ConSeqUMI-\d{8}-\d{6}.*"

def test__conseq__set_command_line_settings(parser): pass

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
    errorOutput = f"argument command: invalid choice: '{nonExistentCommand}' (choose from 'umi')"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args([nonExistentCommand])

def test__conseq__set_command_line_settings__umi_fails_when_does_not_include_input_directory(parser, files, umiArgs):
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
    parentDirectoryPath = "/".join(args["output"].split("/")[:-1])
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
