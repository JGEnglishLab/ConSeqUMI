from PyQt5.QtWidgets import (
    #QLabel,
    #QLineEdit,
    #QWidget,
    #QApplication,
    #QFormLayout,
    #QComboBox,
    #QCheckBox,
    #QPushButton,
    #QFileDialog,
    #QPlainTextEdit,
    #QVBoxLayout,
    #QStyle,
    QMainWindow,
)
from PyQt5.QtCore import Qt
from gui.TabWindow import TabWindow
from gui.TableWindow import TableWindow

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.tableWindow = TableWindow(self)
        self.setCentralWidget(self.tableWindow)
