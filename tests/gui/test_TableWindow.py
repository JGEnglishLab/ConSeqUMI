from PyQt5.QtWidgets import (
    # QLabel,
    # QLineEdit,
    # QWidget,
    QApplication,
    # QFormLayout,
    # QComboBox,
    # QCheckBox,
    # QPushButton,
    # QFileDialog,
    # QPlainTextEdit,
    # QVBoxLayout,
    # QStyle,
    # QMainWindow,
)
import pytest
import sys
import os

srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)

from gui.TableWindow import TableWindow


@pytest.fixture
def tableWindow():
    app = QApplication(sys.argv)
    return TableWindow()


def test__gui_tab_window__initialization(tableWindow):
    assert "tabs" in tableWindow.__dict__
    assert "umiTab" in tableWindow.__dict__
    assert "consTab" in tableWindow.__dict__
    assert "benchmarkTab" in tableWindow.__dict__
