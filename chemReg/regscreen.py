import sys, dbInterface, os, shutil, logging, re
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIntValidator
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QFileDialog, QMessageBox

from chemreglib import *
from sdfreg import LoadSDF
from addmetatags import AddMetaTags


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
        self.addmeta_btn.clicked.connect(self.gotoAddMeta)
        self.regcompound_btn.clicked.connect(self.regCompound)
        self.window().setWindowTitle("Register new compound")
        self.onlyInt = QIntValidator()
        self.populated = False
        self.regcompound_btn.setEnabled(False)
        self.batchOk = False
        self.molOK = False
        #### Regno
        if regno == None:
            self.regno = dbInterface.getNextRegno(self.token)
            dbInterface.createNewRegno(self.regno, self.token)
        else:
            self.regno = str(regno)
        self.regno_eb.setText(str(self.regno))

        self.updateScreen()
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
        self.batch_eb.textChanged.connect(self.batchChanged)
        if regno != None:
            self.batchChanged()
            self.molChanged()
        

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

        self.solvent_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'solvent'))

        self.loadmol_btn.clicked.connect(self.uploadMolfile)
        self.editmol_btn.clicked.connect(self.editMolFile)
        self.editmol_btn.setEnabled(True)

        #self.editmol_btn.setEnabled(False)
        #self.structure_lab.mouseReleaseEvent = self.editMolFile

    def changeLibraryName(self):
        library_name = dbInterface.getLibraryName(self.token,
                                                  self.libraryid_cb.currentText())
        self.librarydesc_eb.setText(library_name)

    def uploadMolfile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                                '.', "Molfiles (*.mol)")
        if fname[0] != '':
            res, sMessage = postMolFile(self, fname[0], self.regno, logging.getLogger(self.mod_name))
            displayMolfile(self)
            updateMoleculeProperties(self)
            self.molChanged()      

    def editMolFile(self, event=None):
        self.editmol_btn.setEnabled(False)
        self.fname = "tmp.mol" # temp file name
        self.fname_path = resource_path(self.fname) # file location / actual file name
        if self.structure_lab.pixmap().isNull():
            # no loaded mol, copying template
            shutil.copy(resource_path("nostruct.mol"), self.fname_path)
        else:
            # download molfile for selected regno, write to 'tmp.mol'
            tmp_mol_str = dbInterface.getMolFile(self.token,
                                                 self.regno)
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
            postMolFile(self, self.fname_path, self.regno, logging.getLogger(self.mod_name))
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

    def allDataPresent(self):

        self.submitter_cb.setProperty("ok", self.submitter_cb.currentText() != ' ')
        self.submitter_cb.style().polish(self.submitter_cb)

        self.compoundtype_cb.setProperty("ok", self.compoundtype_cb.currentText() != ' ')
        self.compoundtype_cb.style().polish(self.compoundtype_cb)

        self.project_cb.setProperty("ok", self.project_cb.currentText() != ' ')
        self.project_cb.style().polish(self.project_cb)

        self.product_cb.setProperty("ok", self.product_cb.currentText() != ' ')
        self.product_cb.style().polish(self.product_cb)

        self.libraryid_cb.setProperty("ok", self.libraryid_cb.currentText() != ' ')
        self.libraryid_cb.style().polish(self.libraryid_cb)

        self.batch_eb.setProperty("ok", self.batchOk)
        self.batch_eb.style().polish(self.batch_eb)

        self.ip_rights_cb.setProperty("ok", self.ip_rights_cb.currentText() != ' ')
        self.ip_rights_cb.style().polish(self.ip_rights_cb)

        self.structure_lab.setProperty("ok", self.molOK)
        self.structure_lab.style().polish(self.structure_lab)
        self.editmol_btn.setProperty("ok", self.molOK)
        self.editmol_btn.style().polish(self.editmol_btn)
        
        if self.submitter_cb.currentText() == ' ' or \
           self.compoundtype_cb.currentText() == ' ' or \
           self.project_cb.currentText() == ' ' or \
           self.product_cb.currentText() == ' ' or \
           self.libraryid_cb.currentText() == ' ' or \
           self.batchOk == False or \
           self.molOK == False or \
           self.ip_rights_cb.currentText() == ' ' or \
           self.compoundid_lab.text() != '':
            return False
        else:
            return True
    
    def batchChanged(self):
        sBatch = self.batch_eb.text()
        pattern = '^[a-zA-Z]{2}[0-9]{7}$'
        self.batchOk = False
        if len(re.findall(pattern, sBatch)) == 1:
            res = dbInterface.updateBatch(self.token,
                                          self.regno,
                                          sBatch)
            if (res == False) and (self.compoundid_lab.text() == ''):
                send_msg("Batch Id error", f"That batch is already in use")
            else:
                self.batchOk = True

            if self.allDataPresent():
                self.regcompound_btn.setEnabled(True)
                return
        self.regcompound_btn.setEnabled(False)
            
    def molChanged(self):
        if (not self.structure_lab.pixmap().isNull()):
            self.molOK = True
        else:
            self.molOK = False
        
        if self.allDataPresent():
            self.regcompound_btn.setEnabled(True)
        else:
            self.regcompound_btn.setEnabled(False)

        
    def changeEvent(self, value='', action='doNothing'):
        if action == 'doNothing':
            return
        self.dirty = True
        if self.allDataPresent():
            self.regcompound_btn.setEnabled(True)
        else:
            self.regcompound_btn.setEnabled(False)
        dbInterface.updateValue(action,
                                value,
                                self.token,
                                self.regno)

    def regCompound(self):
        from searchscreen import SearchScreen
        compound_id = dbInterface.bcpvsRegCompound(self.token,
                                                   self.regno)
        search = SearchScreen(self.token, self.regno)
        self.window().addWidget(search)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

    def gotoSearch(self):
        from searchscreen import SearchScreen
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno,
                                    self.token)
        search = SearchScreen(self.token, self.regno)
        self.window().addWidget(search)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

    def gotoLoadSdf(self):
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno,
                                    self.token)
        loadSDF = LoadSDF(self.token)
        
    def gotoAddMeta(self):
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno,
                                    self.token)
        addMetaTags = AddMetaTags(self.token)

    def updateScreen(self):
        if self.populated == False:
            self.regno_eb.setText(str(self.regno))
            submitters = dbInterface.getColComboData(self.token,
                                                    'chemist')
            self.submitter_cb.addItems(submitters)

            projects = dbInterface.getColComboData(self.token,
                                                'project')
            self.project_cb.addItems(projects)

            compoundTypes = dbInterface.getColComboData(self.token,
                                                        'compound_type')
            self.compoundtype_cb.addItems(compoundTypes)

            productTypes = dbInterface.getColComboData(self.token,
                                                    'product')
            self.product_cb.addItems(productTypes)
            
            libraryIds = dbInterface.getColComboData(self.token,
                                                    'library_id')
            self.libraryid_cb.addItems(libraryIds)
            
            self.ip_rights_cb.addItems(ip_rights_list)

            solvents = dbInterface.getColComboData(self.token,
                                                'solvent')
            self.solvent_cb.addItems(solvents)
        
        # Set current values for this regno
        if self.regno != None:
            self.regno_eb.setText(str(self.regno))
            
            createdDate = dbInterface.getTextColumn(self.token,
                                                    'rdate',
                                                    self.regno)
            self.date_lab.setText(createdDate)

            updateMoleculeProperties(self)
            
            externalCompoundId = dbInterface.getTextColumn(self.token,
                                                        'EXTERNAL_ID',
                                                        self.regno)
            self.externalid_eb.setText(externalCompoundId)

            externalBatch = dbInterface.getTextColumn(self.token,
                                                    'SUPPLIER_BATCH',
                                                    self.regno)
            self.externalbatch_eb.setText(externalBatch)

            currentBatch = dbInterface.getTextColumn(self.token,
                                                    'JPAGE',
                                                    self.regno)
            self.batch_eb.setText(currentBatch)

            compound_id = dbInterface.getTextColumn(self.token,
                                                    'compound_id',
                                                    self.regno)
            self.compoundid_lab.setText(compound_id)

            if compound_id not in (' ', None, ''):
                self.batch_eb.setReadOnly(True)
                self.regcompound_btn.setEnabled(False)
                self.editmol_btn.setEnabled(False)
                self.loadmol_btn.setEnabled(False)
            
            submitter = dbInterface.getTextColumn(self.token,
                                                'chemist',
                                                self.regno)
            self.submitter_cb.setCurrentText(submitter)
            
            project = dbInterface.getTextColumn(self.token,
                                                'project',
                                                self.regno)
            self.project_cb.setCurrentText(project)

            product = dbInterface.getTextColumn(self.token,
                                                'product',
                                                self.regno)
            self.product_cb.setCurrentText(product)

            compoundType = dbInterface.getTextColumn(self.token,
                                                    'compound_type',
                                                    self.regno)
            self.compoundtype_cb.setCurrentText(compoundType)

            libraryId = dbInterface.getTextColumn(self.token,
                                                'library_id',
                                                self.regno)
            self.libraryid_cb.setCurrentText(libraryId)

            library_name = dbInterface.getLibraryName(self.token,
                                                    self.libraryid_cb.currentText())
            self.librarydesc_eb.setText(library_name)

            chromText = dbInterface.getTextColumn(self.token,
                                                'chrom_text',
                                                self.regno)
            self.chrom_text.setPlainText(chromText)
            
            nmrText = dbInterface.getTextColumn(self.token,
                                                'nmr_text',
                                                self.regno)
            self.nmr_text.setPlainText(nmrText)

            commentsText = dbInterface.getTextColumn(self.token,
                                                    'comments',
                                                    self.regno)
            self.comments_text.setPlainText(commentsText)

            msText = dbInterface.getTextColumn(self.token,
                                            'ms_text',
                                            self.regno)
            self.ms_text.setPlainText(msText)

            chromPurity = dbInterface.getTextColumn(self.token,
                                                    'chrom_purity',
                                                    self.regno)
            self.chrompurity_eb.setText(chromPurity)
            nmrPurity = dbInterface.getTextColumn(self.token,
                                                'nmr_purity',
                                                self.regno)
            self.nmrpurity_eb.setText(chromPurity)
            msPurity = dbInterface.getTextColumn(self.token,
                                                'ms_purity',
                                                self.regno)
            self.mspurity_eb.setText(msPurity)
            purity = dbInterface.getTextColumn(self.token,
                                            'purity',
                                            self.regno)
            self.purity_eb.setText(purity)
            
            ip_rights = dbInterface.getTextColumn(self.token,
                                                'ip_rights', 
                                                self.regno)
            self.ip_rights_cb.setCurrentText(ip_rights)

            solvent = dbInterface.getTextColumn(self.token,
                                                'solvent',
                                                self.regno)
            self.solvent_cb.setCurrentText(solvent)

            displayMolfile(self)
        else:
            self.editregno_btn.setEnabled(False)
            self.regno_eb.setText(None)
            self.compoundid_lab.setText(None)
            self.batch_eb.setText(None)
            self.submitter_cb.setCurrentText(' ')
            self.project_cb.setCurrentText(' ')
            self.product_cb.setCurrentText(' ')
            self.compoundtype_cb.setCurrentText(' ')
            self.libraryid_cb.setCurrentText(' ')
            self.chrom_text.setPlainText('')
            self.nmr_text.insertPlainText('')
            self.ms_text.insertPlainText('')
            self.chrompurity_eb.setText(None)
            self.nmrpurity_eb.setText(None)
            self.mspurity_eb.setText(None)
            self.purity_eb.setText(None)
            self.ip_rights_cb.setCurrentText(' ')
            self.externalid_eb.setText(None)
            self.externalbatch_eb.setText(None)
            self.iupac_eb.setText(None)
            self.monoisomass_lab.setText(None)
            self.avgmolmass_lab.setText(None)
            self.date_lab.setText(None)
            self.structure_lab.clear()
            self.solvent_cb.setCurrentText(' ')