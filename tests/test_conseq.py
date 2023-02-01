import pytest
import subprocess
import argparse
import re
from unittest.mock import Mock
from tempfile import TemporaryDirectory, NamedTemporaryFile
import atexit

import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
#testsPath = os.getcwd().split("/")[:-1]
#testsPath = "/".join(testsPath) + "/tests"
#sys.path.insert(1, testsPath)
import conseq

def run_command(commandLine):
    commands = commandLine.split()
    child = subprocess.Popen(commands)
    child.communicate()

@pytest.fixture
def parser():
    def mock_error(message):
        raise ValueError(message)
    argparse.ArgumentParser.error = Mock(side_effect = mock_error)
    parser = conseq.set_command_line_settings()
    return parser

@pytest.fixture
def files():
    class fileObj:
        def __init__(self):
            self.inputDir = TemporaryDirectory(prefix="conseq_input_test_directory_")
            self.outputDir = TemporaryDirectory(prefix="conseq_output_test_directory_")
            self.adapterFile = NamedTemporaryFile(prefix="conseq_adapter_test_file_", suffix=".txt")
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

def test__conseq__generates_original_parser(parser): pass

def test__conseq__fails_when_incorrect_command_provided(parser):
    nonExistentCommand = "nonExistentCommand"
    errorOutput = f"argument command: invalid choice: '{nonExistentCommand}' (choose from 'umi')"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        args = parser.parse_args([nonExistentCommand])

def test__conseq__fails_when_umi_does_not_include_input_directory(parser, files, umiArgs):
    umiArgsWithoutInput = [umiArgs[0]] + umiArgs[3:]
    errorOutput = "the following arguments are required: -i/--input"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithoutInput)

def test__conseq__fails_when_umi_does_not_include_output_directory(parser, files, umiArgs):
    umiArgsWithoutInput = umiArgs[:3] + umiArgs[5:]
    errorOutput = "the following arguments are required: -o/--output"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithoutInput)

def test__conseq__fails_when_umi_does_not_include_adapter_file(parser, files, umiArgs):
    umiArgsWithoutInput = umiArgs[:5]
    errorOutput = "the following arguments are required: -a/--adapter"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        args = parser.parse_args(umiArgsWithoutInput)

    
