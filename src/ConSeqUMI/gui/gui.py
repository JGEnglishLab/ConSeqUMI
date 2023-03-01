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
import sys
import os

from ConSeqUMI.gui.MainWindow import MainWindow

def main():
    app = QApplication(sys.argv)
    view = MainWindow()
    view.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()