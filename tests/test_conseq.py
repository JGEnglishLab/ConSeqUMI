import pytest
import subprocess
import argparse
import re
from unittest.mock import Mock
from tempfile import TemporaryDirectory, NamedTemporaryFile
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
def files(exampleForwardRecord, exampleReverseRecord):
    class fileObj:
        def __init__(self):
            self.inputDir = TemporaryDirectory(prefix="conseq_input_test_directory_")
            self.outputDir = TemporaryDirectory(prefix="conseq_output_test_directory_")
            self.adapterFile = NamedTemporaryFile(prefix="conseq_adapter_test_file_", suffix=".txt")
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
        files.inputDir.name,
        "-o",
        files.outputDir.name,
        "-a",
        files.adapterFile.name,
    ]
    return umiArgs

def test__conseq__set_command_line_settings(parser): pass

def test__conseq__set_command_line_settings__umi_command_succeeds(parser, umiArgs, exampleForwardRecord, exampleReverseRecord):
    args = vars(parser.parse_args(umiArgs))
    assert set([args["input"][0].seq, args["input"][1].seq]) == set([exampleForwardRecord.seq, exampleReverseRecord.seq])
    assert args["output"] == umiArgs[4]
    assert args["adapters"] == umiArgs[6]

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
    umiArgsWithEmptyInput = umiArgs[:2] + [inputFile.name] + umiArgs[3:]
    errorOutput = "The -i or --input argument must be a directory, not a file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithEmptyInput)

def test__conseq__set_command_line_settings__umi_input_directory_fails_when_empty(parser, umiArgs):
    emptyInputDirectory = TemporaryDirectory(prefix="conseq_input_test_directory_empty_")
    umiArgsWithEmptyInput = umiArgs[:2] + [emptyInputDirectory.name] + umiArgs[3:]
    errorOutput = "The -i or --input argument directory must not be empty."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithEmptyInput)

def test__conseq__set_command_line_settings__umi_input_directory_fails_when_contains_non_fastq_file(parser, umiArgs):
    inputDummyFile = NamedTemporaryFile(prefix="conseq_adapter_test_dummy_", dir=umiArgs[2], suffix=".txt", delete=False)
    errorOutput = "The -i or --input argument directory must only contain fastq files (.fq or .fastq)."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgs)
