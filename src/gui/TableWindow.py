from PyQt5.QtWidgets import (
    #QLabel,
    #QLineEdit,
    QWidget,
    #QApplication,
    #QFormLayout,
    #QComboBox,
    #QCheckBox,
    #QPushButton,
    #QFileDialog,
    #QPlainTextEdit,
    QVBoxLayout,
    #QStyle,
    #QMainWindow,
    QTabWidget,
)
from gui.UmiTabWindow import UmiTabWindow
from gui.ConsensusTabWindow import ConsensusTabWindow
from gui.BenchmarkTabWindow import BenchmarkTabWindow

class TableWindow(QWidget):
    def __init__(self, parent=None):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)

        self.tabs = QTabWidget()
        self.umiTab = UmiTabWindow()
        self.consTab = ConsensusTabWindow()
        self.benchmarkTab = BenchmarkTabWindow()
        self.tabs.addTab(self.umiTab,'UMI Processing')
        self.tabs.addTab(self.consTab,'Consensus Generation')
        self.tabs.addTab(self.benchmarkTab,'Benchmarking Consensus')
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

