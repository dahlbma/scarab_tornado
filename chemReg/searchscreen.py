import sys, dbInterface, os, logging
from PyQt5.uic import loadUi
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtCore import QTimer

from chemreglib import *
from sdfreg import LoadSDF
from addmetatags import AddMetaTags

class SearchScreen(QMainWindow):
    def __init__(self, token, regno=None):
        super(SearchScreen, self).__init__()
        self.token = token
        self.mod_name = "search"
        logger = logging.getLogger(self.mod_name)
        loadUi(resource_path("assets/searchchem.ui"), self)
        self.gotoreg_btn.clicked.connect(self.gotoReg)
        self.window().setWindowTitle("Search")
        self.dirty = False
        self.populated = False
        self.searchingInProgress = False
        
        self.regnos = None
        self.editregno_btn.setEnabled(False)
        self.editregno_btn.clicked.connect(self.gotoEditRegno)
        self.openmol_btn.clicked.connect(self.openMolFile)
        
        self.loadsdf_btn.clicked.connect(self.gotoLoadSdf)
        self.back_btn.clicked.connect(self.previousRegno)
        self.forward_btn.clicked.connect(self.nextRegno)
        self.addmeta_btn.clicked.connect(self.gotoAddMeta)
        
        if regno == None:
            self.regno = None
            self.updateScreen()
        else:
            self.searchEvent(regno, action='regno')
        
        self.populated = True

        self.regno_search_eb.textChanged.connect(
            lambda: self.clear_search_fields('regno'))

        self.submitter_search_cb.currentTextChanged.connect(
            lambda: self.clear_search_fields('submitter'))

        self.compoundid_search_eb.textChanged.connect(
            lambda: self.clear_search_fields('compoundid'))

        self.batch_search_eb.textChanged.connect(
            lambda: self.clear_search_fields('batch'))

        self.search_search_btn.clicked.connect(self.launch_searchEvent)

    def openMolFile(self):
        def enable_btn(self):
            self.openmol_btn.setEnabled(True)
        self.openmol_btn.setEnabled(False)
        QTimer.singleShot(5000, lambda: enable_btn(self))
        fname = "tmp.mol" # temp file name
        fname_path = resource_path(fname) # file location / actual file name
        if self.regno != None:
            # download molfile for selected regno, write to 'tmp.mol'
            tmp_mol_str = dbInterface.getMolFile(self.token,
                                                 self.regno)
            tmp_file = open(fname_path, "w")
            n = tmp_file.write(tmp_mol_str)
            tmp_file.close()
            open_file(fname_path)

    

    def gotoLoadSdf(self):
        loadSDF = LoadSDF(self.token)

    def previousRegno(self):
        if self.regno == None or self.regnos == None:
            return
        self.searchingInProgress = True
        newIndex = self.regnos.index(self.regno) -1
        if newIndex > -1:
            self.regno = str(self.regnos[newIndex])
            self.numberOfHits_lab.setText(str(newIndex +1) + " of " + str(len(self.regnos)))
            self.populated == False
            self.updateScreen()
            self.populated == True
        self.searchingInProgress = False
        
    def nextRegno(self):
        if self.regno == None or self.regnos == None:
            return
        self.searchingInProgress = True
        newIndex = self.regnos.index(self.regno) +1
        if newIndex < len(self.regnos):
            self.numberOfHits_lab.setText(str(newIndex +1) + " of " + str(len(self.regnos)))
            self.regno = str(self.regnos[newIndex])
            self.populated == False
            self.updateScreen()
            self.populated == True
        self.searchingInProgress = False

    def gotoEditRegno(self):
        from regscreen import RegScreen
        reg = RegScreen(self.token, self.regno)
        self.window().addWidget(reg)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)
        
    def launch_searchEvent(self):
        #print(f"regno: {self.regno_search_eb.text()}")
        #print(f"submitter: {self.submitter_search_cb.currentText()}")
        #print(f"compoundid: {self.compoundid_search_eb.text()}")
        #print(f"batch: {self.batch_search_eb.text()}")
        if self.regno_search_eb.text() != "":
            self.searchEvent(self.regno_search_eb.text(), 'regno')
        elif self.submitter_search_cb.currentText() != ' ':
            self.searchEvent(self.submitter_search_cb.currentText(), 'CHEMIST')
        elif self.compoundid_search_eb.text() != "":
            self.searchEvent(self.compoundid_search_eb.text(), 'COMPOUND_ID')
        elif self.batch_search_eb.text() != "":
            self.searchEvent(self.batch_search_eb.text(), 'JPAGE')
        else:
            self.searchEvent(None, 'doNothing')

    def searchEvent(self, value='', action='doNothing'):
        if action == 'doNothing':
            self.regnos = []
            self.regno = None
            self.searchStr_lab.setText("No Search")
            self.numberOfHits_lab.setText("0 of 0")
            self.updateScreen()
            return
        # Call search here
        if self.searchingInProgress == False:
            self.searchingInProgress = True
            self.searchCol = action
            self.searchValue = value
            self.searchStr = action + ' = ' + value
            self.searchStr_lab.setText(self.searchStr)
            self.regnos = []
            tmpInts = dbInterface.searchValue(action,
                                              value,
                                              self.token)
            for tmpRegno in tmpInts:
                self.regnos.append(str(tmpRegno))
            if len(self.regnos) > 0 and self.regnos[0] != 'None':
                sString = "1 of " + str(len(self.regnos))
                self.numberOfHits_lab.setText(sString)
                self.regno = str(self.regnos[0])
            else:
                sString = "0 of 0"
                self.numberOfHits_lab.setText(sString)
                self.regno = None
            self.updateScreen()
            self.searchingInProgress = False

    def gotoReg(self):
        from regscreen import RegScreen
        reg = RegScreen(self.token)
        self.window().addWidget(reg)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

    def gotoAddMeta(self):
        addMetaTags = AddMetaTags(self.token)

    def clear_search_fields(self, field):
        self.regno_search_eb.textChanged.disconnect()
        self.submitter_search_cb.currentTextChanged.disconnect()
        self.compoundid_search_eb.textChanged.disconnect()
        self.batch_search_eb.textChanged.disconnect()
        if field == 'regno':
            #self.regno_search_eb.setText(None)
            self.submitter_search_cb.setCurrentText(' ')
            self.compoundid_search_eb.setText(None)
            self.batch_search_eb.setText(None)
        elif field == 'submitter':
            self.regno_search_eb.setText(None)
            #self.submitter_search_cb.setCurrentText(' ')
            self.compoundid_search_eb.setText(None)
            self.batch_search_eb.setText(None)
        elif field == 'compoundid':
            self.regno_search_eb.setText(None)
            self.submitter_search_cb.setCurrentText(' ')
            #self.compoundid_search_eb.setText(None)
            self.batch_search_eb.setText(None)
        elif field == 'batch':
            self.regno_search_eb.setText(None)
            self.submitter_search_cb.setCurrentText(' ')
            self.compoundid_search_eb.setText(None)
            #self.batch_search_eb.setText(None)
        else:
            self.regno_search_eb.setText(None)
            self.submitter_search_cb.setCurrentText(' ')
            self.compoundid_search_eb.setText(None)
            self.batch_search_eb.setText(None)
        self.regno_search_eb.textChanged.connect(
            lambda: self.clear_search_fields('regno'))
        self.submitter_search_cb.currentTextChanged.connect(
            lambda: self.clear_search_fields('submitter'))
        self.compoundid_search_eb.textChanged.connect(
            lambda: self.clear_search_fields('compoundid'))
        self.batch_search_eb.textChanged.connect(
            lambda: self.clear_search_fields('batch'))

    def updateScreen(self):
        if self.populated == False:
            self.regno_eb.setText(str(self.regno))
            submitters = dbInterface.getColComboData(self.token,
                                                     'chemist')
            self.submitter_cb.addItems(submitters)
            try:
                self.submitter_search_cb.clear()
                self.submitter_search_cb.addItems(submitters)
            except:
                pass
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
        if self.regno not in (None, 'None'):
            self.editregno_btn.setEnabled(True)
            self.openmol_btn.setEnabled(True)
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
            if libraryId == None:
                libraryId = ' '
            self.libraryid_cb.setCurrentText(libraryId)

            library_name = dbInterface.getLibraryName(self.token,
                                                      libraryId)
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
            self.submitter_cb.setCurrentText(None)
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
            self.chns_eb.setText(None)
            self.mf_eb.setText(None)
            self.librarydesc_eb.setText(None)
            self.date_lab.setText(None)
            self.structure_lab.clear()
            self.solvent_cb.setCurrentText(' ')
            self.openmol_btn.setEnabled(False)
