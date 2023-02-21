from PyQt5.QtWidgets import (
    #QLabel,
    #QLineEdit,
    #QWidget,
    QApplication,
    #QFormLayout,
    #QComboBox,
    #QCheckBox,
    #QPushButton,
    #QFileDialog,
    #QPlainTextEdit,
    #QVBoxLayout,
    #QStyle,
    #QMainWindow,
)
import pytest
import sys
import os
srcPath = os.getcwd().split("/")[:-1]
srcPath = "/".join(srcPath) + "/src"
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
    assert "tab1" in tableWindow.__dict__
    assert "tab2" in tableWindow.__dict__
    assert "tab3" in tableWindow.__dict__

