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
        self.dirty = False
        self.updateStructure_btn.clicked.connect(self.updateStructure)
        self.window().setWindowTitle("Edit structure")
        self.updateStructure_btn.setEnabled(False)
        self.batchOk = False
        self.molOK = False
        self.compound_id = ''
        
        self.compoundId_eb.textChanged.connect(self.compoundIdChanged)

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
        self.fs_watcher.fileChanged.connect(self.test_msg)
        
        return

        
        if self.new_structure_lab.pixmap().isNull():
            # no loaded mol, copying template
            shutil.copy(resource_path("nostruct.mol"), self.fname_path)
        else:
            # download molfile for selected regno, write to 'tmp.mol'
            tmp_mol_str = dbInterface.getMolFileBcpvs(self.token,
                                                 self.compound_id)
            tmp_file = open(self.fname_path, "w")
            n = tmp_file.write(tmp_mol_str)
            tmp_file.close()
        retcode = open_file(self.fname_path)
        # confirm dialogue
        self.fs_watcher = QtCore.QFileSystemWatcher([self.fname_path])
        self.fs_watcher.fileChanged.connect(self.test_msg)

    def test_msg(self):
        self.fs_watcher = None
        msg = QMessageBox()
        msg.setWindowTitle("Edit " + self.fname)
        msg.setIcon(QMessageBox.Question)
        msg.setText("Do you want to save .mol-file changes to the database?")
        msg.setStandardButtons(QMessageBox.Save | QMessageBox.Cancel)
        msg.setDefaultButton(QMessageBox.Save)
        btnS = msg.button(QMessageBox.Save)
        msg.exec_()
        ok_msg = QMessageBox()
        ok_msg.setStandardButtons(QMessageBox.Ok)
        ok_msg.setWindowTitle("Edit " + self.fname)
        if (msg.clickedButton() == btnS):
            # save changes
            #postMolFile(self, self.fname_path, self.regno, logging.getLogger(self.mod_name))
            displayMolfile(self)
            updateMoleculeProperties(self)
            ok_msg.setText("Updated .mol file")
            self.molChanged()
        else:
            # cancel, do nothing
            ok_msg.setText("Did not update .mol file in database. 'tmp.mol' deleted.")
        ok_msg.exec_()
        # cleanup, remove tmp.mol
        os.remove(self.fname_path)
        self.fname = None
        self.fname_path = None
        self.editmol_btn.setEnabled(True)

    def compoundIdChanged(self):
        sCompoundId = self.compoundId_eb.text()
        pattern = '^CBK[0-9]{6}$'
        self.batchOk = False
        if len(re.findall(pattern, sCompoundId)) == 1:
            self.compound_id = sCompoundId
            res = dbInterface.getMolFileBcpvs(self.token, sCompoundId)
        self.updateStructure_btn.setEnabled(False)
        
    def molChanged(self):
        if (not self.new_structure_lab.pixmap().isNull()):
            self.molOK = True
        else:
            self.molOK = False
            self.molOK = True
        
        if self.allDataPresent():
            self.regcompound_btn.setEnabled(True)
        else:
            self.regcompound_btn.setEnabled(False)

    def updateStructure(self):
        pass
        #compound_id = dbInterface.bcpvsRegCompound(self.token,
        #                                           self.regno)
        #search = SearchScreen(self.token, self.regno)
        #self.window().addWidget(search)
        #self.window().setCurrentIndex(self.window().currentIndex() + 1)
