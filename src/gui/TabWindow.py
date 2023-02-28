from PyQt5.QtWidgets import (
    QLabel,
    #QLineEdit,
    QWidget,
    #QApplication,
    QFormLayout,
    #QComboBox,
    #QCheckBox,
    QPushButton,
    #QFileDialog,
    QPlainTextEdit,
    QVBoxLayout,
    #QStyle,
    #QMainWindow,
)
from abc import ABC, ABCMeta, abstractmethod

class AbstractTabWindow(ABCMeta, type(QWidget)):
    pass

class TabWindow(QWidget, metaclass=AbstractTabWindow):
    def __init__(self):
        super().__init__()

        self.p = None

        dlgLayout = QVBoxLayout()
        fileLayout = QFormLayout()
        settingLayout = QFormLayout()

        self.set_file_layout(fileLayout)
        dlgLayout.addWidget(QLabel('File Input:'))
        dlgLayout.addLayout(fileLayout)

        self.set_setting_layout(settingLayout)
        dlgLayout.addWidget(QLabel('Settings:'))
        dlgLayout.addLayout(settingLayout)

        self.runBtn = QPushButton('Execute')
        self.killBtn = QPushButton('Kill Process')
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        dlgLayout.addWidget(self.runBtn)
        dlgLayout.addWidget(self.killBtn)
        dlgLayout.addWidget(self.text)
        self.runBtn.setDefault(True)
        self.killBtn.setEnabled(False)
        self.killBtn.clicked.connect(self.kill_process)
        #self.runBtn.clicked.connect(self.start_process)

        self.setLayout(dlgLayout)

    @abstractmethod
    def set_file_layout(self, fileLayout: type(QFormLayout)) -> None:
        pass

    @abstractmethod
    def set_setting_layout(self, settingLayout: QFormLayout) -> None:
        pass

    def message(self, s):
        self.text.appendPlainText(s)

    def kill_process(self):
        self.p.kill()

    def process_finished(self):
        self.message("Process finished.")
        self.killBtn.setEnabled(False)
        self.p = None