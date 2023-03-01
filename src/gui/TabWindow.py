from PyQt5.QtWidgets import (
    QLabel,
    #QLineEdit,
    QWidget,
    #QApplication,
    QFormLayout,
    #QComboBox,
    #QCheckBox,
    QPushButton,
    QFileDialog,
    QPlainTextEdit,
    QVBoxLayout,
    #QStyle,
    #QMainWindow,
)
from PyQt5.QtCore import QProcess
from PyQt5.QtGui import QFont
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
        fileInputLabel = QLabel('File Input:')
        fileInputLabel.setFont(QFont("Times", 17, QFont.Bold))
        dlgLayout.addWidget(fileInputLabel)
        dlgLayout.addLayout(fileLayout)

        self.set_setting_layout(settingLayout)
        settingsInputLabel = QLabel('Settings:')
        settingsInputLabel.setFont(QFont("Times", 17, QFont.Bold))
        dlgLayout.addWidget(settingsInputLabel)
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
        self.runBtn.clicked.connect(self.start_process)

        self.setLayout(dlgLayout)

    @abstractmethod
    def set_file_layout(self, fileLayout: type(QFormLayout)) -> None:
        pass

    @abstractmethod
    def set_setting_layout(self, settingLayout: QFormLayout) -> None:
        pass

    @abstractmethod
    def set_args(self) -> list:
        pass

    def get_file(self, text, isFile=False):
        dialog = QFileDialog()
        if isFile:
            fname = QFileDialog.getOpenFileName(self, "Open File", "c:\\")
            text.setText(fname[0])
        else:
            fname = QFileDialog.getExistingDirectory(self, "Open Directory", "c:\\")
            text.setText(fname)

    def message(self, s):
        self.text.appendPlainText(s)

    def kill_process(self):
        self.p.kill()

    def process_finished(self):
        self.message("Process finished.")
        self.killBtn.setEnabled(False)
        self.p = None

    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)

    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)

    def start_process(self):
        args = self.set_args()
        self.killBtn.setEnabled(True)

        if self.p is None:  # No process running.

            self.message("Executing process")
            self.message("python " + " ".join(args))
            self.p = (
                QProcess()
            )  # Keep a reference to the QProcess (e.g. on self) while it's running.
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.finished.connect(self.process_finished)  # Clean up once complete.
            self.p.start("python", args)