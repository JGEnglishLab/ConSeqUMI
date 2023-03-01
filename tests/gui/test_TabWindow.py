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
srcPath = "/".join(srcPath) + "/src/ConSeqUMI"
sys.path.insert(1, srcPath)
testsPath = os.getcwd().split("/")[:-1]
testsPath = "/".join(testsPath) + "/tests"
sys.path.insert(1, testsPath)

