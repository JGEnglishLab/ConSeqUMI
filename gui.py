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
)

class longreadWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.p = None
        self.dict = None
        self.INVALID = -1
        self.VALID_COLOR = 'lightgreen'
        self.INVALID_COLOR = 'rgb(255,99,71)'

        self.setWindowTitle('longread_umi_python')
        formLayout = QFormLayout()

        self.inputLabel = QLabel('Input Directory Path')
        self.inputField = QLineEdit()
        self.inputField.setEnabled(False)
        self.inputBrowseButton = QPushButton('Browse')
        self.inputBrowseButton.clicked.connect(lambda:self.get_file(self.inputField))

        formLayout.addRow(self.inputLabel, self.inputBrowseButton)
        formLayout.addRow(self.inputField)


        self.outputLabel = QLabel('Output Directory Path')
        self.outputField = QLineEdit()
        self.outputField.setEnabled(False)
        self.outputBrowseButton = QPushButton('Browse')
        self.outputBrowseButton.clicked.connect(lambda:self.get_file(self.outputField))

        formLayout.addRow(self.outputLabel, self.outputBrowseButton)
        formLayout.addRow(self.outputField)

        self.outputNameLabel = QLabel('Output Directory Title')
        self.outputNameField = QLineEdit()
        formLayout.addRow(self.outputNameLabel, self.outputNameField)


        self.adapterLabel = QLabel('Adapter File Path')
        self.adapterField = QLineEdit()
        self.adapterField.setEnabled(False)
        self.adapterBrowseButton = QPushButton('Browse')
        self.adapterBrowseButton.clicked.connect(lambda:self.get_file(self.adapterField, isFile=True))

        formLayout.addRow(self.adapterLabel, self.adapterBrowseButton)
        formLayout.addRow(self.adapterField)

        self.variantLabel = QLabel('Consolidate Variants in Final Output?')
        self.variantCheckBox = QCheckBox()
        formLayout.addRow(self.variantLabel, self.variantCheckBox)

        self.runBtn = QPushButton('Execute')
        self.killBtn = QPushButton('Kill Process')
        self.processText = QPlainTextEdit()
        self.processText.setReadOnly(True)
        formLayout.addWidget(self.runBtn)
        formLayout.addWidget(self.killBtn)
        formLayout.addWidget(self.processText)
        self.runBtn.setDefault(True)
        self.killBtn.setEnabled(False)


        self.setLayout(formLayout)

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
        tempDict['adapterFile'] = self.return_file_path_value(self.adapterField.text(), self.adapterLabel, permittedTypes=['txt'])
        if self.variantCheckBox.isChecked(): tempDict['isVariant'] = True; self.set_text_color(self.variantLabel, isValid=True)
        else: tempDict['isVariant'] = False; self.set_text_color(self.variantLabel, isValid=False)
        if False in list(tempDict.values()): return False
        args = []
        args += ['i', tempDict['inputDir']]
        args += ['o', tempDict['outputDir']]
        args += ['a', tempDict['adapterFile']]
        if tempDict['isVariant']: args += ['v']
        return args

    def return_file_path_value(self, path, text, permittedTypes=[]):
        if len(path)==0: return self.set_text_color(text)
        if len(permittedTypes)!=0 and path.split('.')[-1].lower() not in permittedTypes: return self.set_text_color(text)
        return self.set_text_color(text, path)

    def start_process(self):
        args = self.set_args()
        if not args: return
        self.killBtn.setEnabled(True)

        if self.p is None:  # No process running.
            args = []

            self.message("Executing process")
            #self.p = QProcess()  # Keep a reference to the QProcess (e.g. on self) while it's running.
            #self.p.readyReadStandardOutput.connect(self.handle_stdout)
            #self.p.readyReadStandardError.connect(self.handle_stderr)
            #self.p.finished.connect(self.process_finished)  # Clean up once complete.
            #self.p.start('csodiaq', args)

application = QApplication([])
mainWindow = longreadWindow()


mainWindow.show()

application.exec()
