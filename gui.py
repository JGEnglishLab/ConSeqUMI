from PyQt5.QtWidgets import (
    QLabel,
    QLineEdit,
    QWidget,
    QApplication,
    QFormLayout,
    QCheckBox,
    QPushButton,
    QFileDialog,
    QPlainTextEdit,
    QVBoxLayout,
)

from PyQt5.QtCore import QProcess

class longreadWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.p = None
        self.dict = None
        self.INVALID = -1
        self.VALID_COLOR = 'lightgreen'
        self.INVALID_COLOR = 'rgb(255,99,71)'

        self.setWindowTitle('longread_umi_python')
        dlgLayout = QVBoxLayout()
        formLayout = QFormLayout()

        self.inputLabel = QLabel('Input Directory Path')
        self.inputField = QLineEdit()
        self.inputField.setEnabled(False)
        self.inputField.setText('/Users/calebcranney/Documents/Projects/longread_umi_python/test/data/7_UMI_test_input/barcode01')
        self.inputBrowseButton = QPushButton('Browse')
        self.inputBrowseButton.clicked.connect(lambda:self.get_file(self.inputField))
        #self.inputCopyButton = QPushButton('Copy to Clipboard')
        #self.inputCopyButton.clicked.connect(lambda:self.copy_to_clipboard(self.inputField.text()))

        formLayout.addRow(self.inputLabel, self.inputBrowseButton)
        formLayout.addRow(self.inputField)


        self.outputLabel = QLabel('Output Directory Path')
        self.outputField = QLineEdit()
        self.outputField.setEnabled(False)
        self.outputField.setText('/Users/calebcranney/Documents/Projects/longread_umi_python/test/data')
        self.outputBrowseButton = QPushButton('Browse')
        self.outputBrowseButton.clicked.connect(lambda:self.get_file(self.outputField))
        #self.outputCopyButton = QPushButton('Copy to Clipboard')
        #self.outputCopyButton.clicked.connect(lambda:self.copy_to_clipboard(self.outputField.text()))


        formLayout.addRow(self.outputLabel, self.outputBrowseButton)
        formLayout.addRow(self.outputField)

        self.outputNameLabel = QLabel('Output Directory Title')
        self.outputNameField = QLineEdit()
        self.outputNameField.setText('delete')
        formLayout.addRow(self.outputNameLabel, self.outputNameField)


        self.adapterLabel = QLabel('Adapter File Path')
        self.adapterField = QLineEdit()
        self.adapterField.setEnabled(False)
        self.adapterField.setText('/Users/calebcranney/Documents/Projects/longread_umi_python/test/data/adapters.txt')
        self.adapterBrowseButton = QPushButton('Browse')
        self.adapterBrowseButton.clicked.connect(lambda:self.get_file(self.adapterField, isFile=True))
        #self.adapterCopyButton = QPushButton('Copy to Clipboard')
        #self.adapterCopyButton.clicked.connect(lambda:self.copy_to_clipboard(self.adapterField.text()))

        formLayout.addRow(self.adapterLabel, self.adapterBrowseButton)
        formLayout.addRow(self.adapterField)

        self.variantLabel = QLabel('Consolidate Variants in Final Output?')
        self.variantCheckBox = QCheckBox()
        formLayout.addRow(self.variantLabel, self.variantCheckBox)


        self.benchmarkLabel = QLabel('Perform Benchmarking Run (note: runtime significantly longer)')
        self.benchmarkCheckBox = QCheckBox()
        formLayout.addRow(self.benchmarkLabel, self.benchmarkCheckBox)
        self.benchmarkCheckBox.stateChanged.connect(lambda:self.disable_second_box(self.benchmarkCheckBox, self.variantCheckBox))

        dlgLayout.addLayout(formLayout)

        self.runBtn = QPushButton('Execute')
        self.killBtn = QPushButton('Kill Process')
        self.processText = QPlainTextEdit()
        self.processText.setReadOnly(True)
        self.processCopyButton = QPushButton('Copy to Clipboard')
        self.processCopyButton.clicked.connect(lambda:self.copy_to_clipboard(self.processText.toPlainText()))

        dlgLayout.addWidget(self.runBtn)
        dlgLayout.addWidget(self.killBtn)
        dlgLayout.addWidget(self.processText)
        dlgLayout.addWidget(self.processCopyButton)
        self.runBtn.setDefault(True)
        self.killBtn.setEnabled(False)
        self.killBtn.clicked.connect(self.kill_process)
        self.runBtn.clicked.connect(self.start_process)

        self.setLayout(dlgLayout)

    def disable_second_box(self, checkBox1, checkBox2):
        if checkBox1.isChecked():
            checkBox2.setCheckState(False)
            checkBox2.setEnabled(False)
        else:
            checkBox2.setEnabled(True)

    def get_file(self, text, isFile=False):
        dialog = QFileDialog()
        if isFile: fname = QFileDialog.getOpenFileName(self, 'Open File', 'c:\\'); text.setText(fname[0])
        else: fname = QFileDialog.getExistingDirectory(self, 'Open Directory', 'c:\\'); text.setText(fname)

    def set_text_color(self, text, isValid=False):
        if isValid: text.setStyleSheet('background-color: ' + self.VALID_COLOR)
        else: text.setStyleSheet('background-color: ' + self.INVALID_COLOR)

    def set_args(self):
        tempDict = {}
        tempDict['inputDir'] = self.return_file_path_value(self.inputField.text(), self.inputLabel)
        tempDict['outputDir'] = self.return_file_path_value(self.outputField.text(), self.outputLabel)
        if len(self.outputNameField.text())==0: tempDict['outputName'] = False; self.set_text_color(self.outputNameLabel)
        else: tempDict['outputName'] = self.outputNameField.text(); self.set_text_color(self.outputNameLabel, isValid = True)
        tempDict['adapterFile'] = self.return_file_path_value(self.adapterField.text(), self.adapterLabel, permittedTypes=['txt'])
        if self.variantCheckBox.isChecked(): tempDict['isVariant'] = 2; self.set_text_color(self.variantLabel, isValid=True)
        else: tempDict['isVariant'] = 1; self.set_text_color(self.variantLabel, isValid=True)
        if self.benchmarkCheckBox.isChecked(): tempDict['isBenchmark'] = 2; self.set_text_color(self.benchmarkLabel, isValid=True)
        else: tempDict['isBenchmark'] = 1; self.set_text_color(self.benchmarkLabel, isValid=True)

        if False in list(tempDict.values()): return False
        args = []
        args += ['pipeline.py']
        args += ['i', tempDict['inputDir']]
        args += ['o', tempDict['outputDir'] + '/' + tempDict['outputName']]
        args += ['a', tempDict['adapterFile']]
        if tempDict['isVariant']==2: args += ['v']
        if tempDict['isBenchmark']==2: args += ['bc']
        return args

    def return_file_path_value(self, path, text, permittedTypes=[]):
        if len(path)==0: self.set_text_color(text); return False
        if len(permittedTypes)!=0 and path.split('.')[-1].lower() not in permittedTypes: self.set_text_color(text); return False
        self.set_text_color(text, isValid=True)
        return path

    def message(self, s):
        self.processText.appendPlainText(s)

    def start_process(self):
        args = self.set_args()
        if not args: return
        self.killBtn.setEnabled(True)

        if self.p is None:  # No process running.

            self.message("Executing process")
            self.message("Input Directory: " + self.inputField.text())
            self.message("Output Directory: " + self.outputField.text() + '/' + self.outputNameField.text())
            self.message("Adapter File: " + self.adapterField.text())
            if self.variantCheckBox.isChecked(): self.message('Variant setting enabled')
            if self.benchmarkCheckBox.isChecked(): self.message('Benchmark setting enabled')
            self.message("\n\n")
            self.p = QProcess()  # Keep a reference to the QProcess (e.g. on self) while it's running.
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.finished.connect(self.process_finished)  # Clean up once complete.
            self.p.start('python', args)

    def kill_process(self):
        self.p.kill()

    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)

    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)

    def process_finished(self):
        self.message("Process finished.")
        self.killBtn.setEnabled(False)
        self.p = None

    def copy_to_clipboard(self, text):
        cb = QApplication.clipboard()
        cb.clear(mode=cb.Clipboard)
        cb.setText(text, mode=cb.Clipboard)
        self.message("content is copied to clipboard")

application = QApplication([])
mainWindow = longreadWindow()


mainWindow.show()

application.exec()
