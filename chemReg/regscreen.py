import sys, dbInterface, os, shutil, logging
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIntValidator
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QFileDialog, QMessageBox

from chemreglib import *
from sdfreg import LoadSDF


class RegScreen(QMainWindow):
    def __init__(self, token, regno=None):
        super(RegScreen, self).__init__()
        loadUi(resource_path("assets/regchem.ui"), self)
        self.token = token
        self.mod_name = "reg"
        logger = logging.getLogger(self.mod_name)
        self.dirty = False
        self.loadsdf_btn.clicked.connect(self.gotoLoadSdf)
        self.gotosearch_btn.clicked.connect(self.gotoSearch)
        self.window().setWindowTitle("Register new compound")
        self.onlyInt = QIntValidator()
        self.populated = False

        #### Regno
        if regno == None:
            self.regno = dbInterface.getNextRegno(self.token)
            dbInterface.createNewRegno(self.regno, self.token)
        else:
            self.regno = regno
        self.regno_eb.setText(self.regno)

        updateScreen(self)
        self.populated = True
        #### Chemist
        self.submitter_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'CHEMIST'))

        #### Projects
        self.project_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'project'))

        #### CompoundTypes
        self.compoundtype_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'compound_type'))

        #### ProductTypes
        self.product_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'product'))

        #### Libraries and LibraryName
        self.libraryid_cb.currentTextChanged.connect(self.changeLibraryName)
        self.libraryid_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'library_id'))

        #### ExternalId
        self.externalid_eb.textChanged.connect(
            lambda: self.changeEvent(self.externalid_eb.text(), 'EXTERNAL_ID'))

        #### ExternalBatch
        self.externalbatch_eb.textChanged.connect(
            lambda: self.changeEvent(self.externalbatch_eb.text(), 'SUPPLIER_BATCH'))

        #### Batch
        self.batch_eb.editingFinished.connect(
            lambda: self.changeEvent(self.batch_eb.text(), 'JPAGE'))

        self.chrom_text.textChanged.connect(
            lambda: self.changeEvent(self.chrom_text.toPlainText(), 'CHROM_TEXT'))
        self.chrompurity_eb.textChanged.connect(
            lambda: self.changeEvent(self.chrompurity_eb.text(), 'CHROM_PURITY'))
        self.chrompurity_eb.setValidator(self.onlyInt)
        
        self.nmr_text.textChanged.connect(
            lambda: self.changeEvent(self.nmr_text.toPlainText(), 'NMR_TEXT'))
        self.nmrpurity_eb.textChanged.connect(
            lambda: self.changeEvent(self.nmrpurity_eb.text(), 'NMR_PURITY'))
        self.nmrpurity_eb.setValidator(self.onlyInt)
        
        self.ms_text.textChanged.connect(
            lambda: self.changeEvent(self.ms_text.toPlainText(), 'MS_TEXT'))
        self.mspurity_eb.textChanged.connect(
            lambda: self.changeEvent(self.mspurity_eb.text(), 'MS_PURITY'))
        self.mspurity_eb.setValidator(self.onlyInt)
        
        self.ip_rights_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'ip_rights'))
        
        self.purity_eb.textChanged.connect(
            lambda: self.changeEvent(self.purity_eb.text(), 'PURITY'))
        self.purity_eb.setValidator(self.onlyInt)
        
        self.comments_text.textChanged.connect(
            lambda: self.changeEvent(self.comments_text.toPlainText(), 'COMMENTS'))

        self.loadmol_btn.clicked.connect(self.uploadMolfile)
        self.editmol_btn.clicked.connect(self.editMolFile)
        self.structure_lab.mouseReleaseEvent = self.editMolFile

    def changeLibraryName(self):
        library_name = dbInterface.getLibraryName(self.token,
                                                  self.libraryid_cb.currentText())
        self.librarydesc_eb.setText(library_name)

    def uploadMolfile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                                '.', "Molfiles (*.mol)")
        if fname[0] != '':
            postMolFile(self, fname[0], self.regno, logging.getLogger(self.mod_name))
            displayMolfile(self)
            updateMoleculeProperties(self)
            
    def editMolFile(self, event=None):
        fname = "tmp.mol" # temp file name
        fname_path = resource_path(fname) # file location / actual file name
        if self.structure_lab.pixmap().isNull():
            # no loaded mol, copying template
            shutil.copy(resource_path("nostruct.mol"), fname_path)
        else:
            # download molfile for selected regno, write to 'tmp.mol'
            tmp_mol_str = dbInterface.getMolFile(self.token, self.regno)
            tmp_file = open(fname_path, "w")
            n = tmp_file.write(tmp_mol_str)
            tmp_file.close()
        retcode = open_file(fname_path)
        # confirm dialogue
        msg = QMessageBox()
        msg.setWindowTitle("Edit " + fname)
        msg.setIcon(QMessageBox.Question)
        msg.setText("Do you want to save .mol-file changes to the database?")
        msg.setStandardButtons(QMessageBox.Save | QMessageBox.Cancel)
        msg.setDefaultButton(QMessageBox.Save)
        btnS = msg.button(QMessageBox.Save)
        msg.exec_()
        ok_msg = QMessageBox()
        ok_msg.setStandardButtons(QMessageBox.Ok)
        ok_msg.setWindowTitle("Edit " + fname)
        if (msg.clickedButton() == btnS):
            # save changes
            postMolFile(self, fname_path, self.regno, logging.getLogger(self.mod_name))
            displayMolfile(self)
            ok_msg.setText("Updated .mol file")
        else:
            # cancel, do nothing
            ok_msg.setText("Did not update .mol file")
        ok_msg.exec_()
        # cleanup, remove tmp.mol
        os.remove(fname_path) 
    
    def changeEvent(self, value='', action='doNothing'):
        if action == 'doNothing':
            return
        self.dirty = True
        dbInterface.updateValue(action, value, self.token, self.regno)

    def gotoSearch(self):
        from searchscreen import SearchScreen
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno, self.token)
        search = SearchScreen(self.token)
        self.window().addWidget(search)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

    def gotoLoadSdf(self):
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno, self.token)
        loadSDF = LoadSDF(self.token)
