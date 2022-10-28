import sys, dbInterface, os, shutil, logging, re
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIntValidator
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QFileDialog, QMessageBox

from chemreglib import *


class EditStructure(QMainWindow):
    def __init__(self, token):
        super(EditStructure, self).__init__()
        loadUi(resource_path("assets/structure.ui"), self)
        self.token = token
        self.mod_name = "reg"
        logger = logging.getLogger(self.mod_name)
        self.updateStructure_btn.clicked.connect(self.updateStructure)
        self.window().setWindowTitle("Edit structure")
        self.updateStructure_btn.setEnabled(False)
        self.compound_id = ''
        self.regno = ''
        
        self.compoundId_eb.textChanged.connect(self.compoundIdChanged)
        self.regno_eb.textChanged.connect(self.regnoChanged)

        self.editStructure_btn.clicked.connect(self.editMolFile)
        self.editStructure_btn.setEnabled(True)


    def uploadMolfile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                                '.', "Molfiles (*.mol)")
        if fname[0] != '':
            #res, sMessage = postMolFile(self, fname[0], self.regno, logging.getLogger(self.mod_name))
            displayMolfile(self)
            updateMoleculeProperties(self)
            self.molChanged()      

    def editMolFile(self, event=None):
        self.editStructure_btn.setEnabled(False)
        self.fname = "tmp.mol" # temp file name
        self.fname_path = resource_path(self.fname) # file location / actual file name
        shutil.copy(resource_path("nostruct.mol"), self.fname_path)
        retcode = open_file(self.fname_path)

        self.fs_watcher = QtCore.QFileSystemWatcher([self.fname_path])
        self.fs_watcher.fileChanged.connect(self.createMolImage)
        

    def createMolImage(self):
        molfile_file = open(self.fname_path, "r")
        molfile = molfile_file.read()
        molfile_file.close()

        dbInterface.createMolImageFromMolfile(self.token, molfile)
        displayMolfile(self, 'new_structure', self.new_structure_lab)
    

    def compoundIdChanged(self):
        sCompoundId = self.compoundId_eb.text().upper()
        pattern = '^CBK[0-9]{6}$'
        if len(re.findall(pattern, sCompoundId)) == 1:
            self.compound_id = sCompoundId
            res = dbInterface.createBcpvsMolImage(self.token, sCompoundId)
            displayMolfile(self, sCompoundId, self.current_structure_lab)
        else:
            self.updateStructure_btn.setEnabled(False)


    def regnoChanged(self):
        sRegno = self.regno_eb.text()
        #3913654
        pattern = '^[0-9]{7}$'
        if len(re.findall(pattern, sRegno)) == 1:
            self.regno = sRegno
            res = dbInterface.createMolImage(self.token, sRegno)
            displayMolfile(self, sRegno, self.original_structure_lab)
        else:
            self.updateStructure_btn.setEnabled(False)


    def molChanged(self):
        if (not self.new_structure_lab.pixmap().isNull()):
            self.molOK = True
        else:
            self.molOK = False
            self.molOK = True
        
    def updateStructure(self):
        pass
        #compound_id = dbInterface.bcpvsRegCompound(self.token,
        #                                           self.regno)
        #search = SearchScreen(self.token, self.regno)
        #self.window().addWidget(search)
        #self.window().setCurrentIndex(self.window().currentIndex() + 1)
