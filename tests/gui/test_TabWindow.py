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

from gui.TabWindow import TabWindow

@pytest.fixture
def tabWindow():
    app = QApplication(sys.argv)
    return TabWindow()

def test__gui_tab_window__initialization(tabWindow): pass
