import pytest
import subprocess
import argparse
import re
from unittest.mock import Mock

import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
import conseq

def run_command(commandLine):
    commands = commandLine.split()
    child = subprocess.Popen(commands)
    child.communicate()

@pytest.fixture
def parser():
    parser = conseq.set_command_line_settings()
    def mock_error(message):
        raise ValueError(message)
    parser.error = Mock(side_effect = mock_error)
    return parser

def test__conseq__generates_original_parser(parser): pass

def test__conseq__fails_when_incorrect_command_provided(parser):
    nonExistentCommand = "nonExistentCommand"
    errorOutput = f"argument command: invalid choice: '{nonExistentCommand}' (choose from 'umi')"
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        args = parser.parse_args([nonExistentCommand])
