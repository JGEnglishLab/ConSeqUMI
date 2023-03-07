import pytest
from unittest.mock import Mock
import re
import sys
import os

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)
from Printer import Printer


@pytest.fixture
def printer():
    printer = Printer()
    printer.startTime = 0.0
    return printer


def test__Printer__initialization():
    printer = Printer()
    assert printer.startTime
    assert type(printer.startTime) == type(0.0)


def test__Printer__determine_time_since_start_in_human_readable_format(printer):
    assert printer.determine_time_since_start_in_human_readable_format(0) == "00:00:00"
    assert printer.determine_time_since_start_in_human_readable_format(1) == "00:00:01"
    assert printer.determine_time_since_start_in_human_readable_format(60) == "00:01:00"
    assert (
        printer.determine_time_since_start_in_human_readable_format(3600) == "01:00:00"
    )


def test__Printer__succeeds_when_called_directly(printer, capsys):
    text = "random text"
    timePattern = r"\d\d:\d\d:\d\d"
    matchPattern = "\n----> " + timePattern + ": " + text + "\n\n"
    printer(text)
    captured = capsys.readouterr()
    assert re.match(matchPattern, captured.out)
