import sys, dbInterface, os, logging
from PyQt5.uic import loadUi
from PyQt5.QtWidgets import QMainWindow

from chemreglib import *

from sdfreg import LoadSDF
from addmetatags import AddMetaTags

class SearchScreen(QMainWindow):
    def __init__(self, token, database, regno=None):
        super(SearchScreen, self).__init__()
        self.token = token
        self.database = database
        self.mod_name = "search"
        logger = logging.getLogger(self.mod_name)
        loadUi(resource_path("assets/searchchem.ui"), self)
        self.gotoreg_btn.clicked.connect(self.gotoReg)
        self.window().setWindowTitle("Search")
        self.dirty = False
        self.populated = False
        self.searchingInProgress = False
        if regno == None:
            self.regno = None
        else:
            self.searchEvent(regno, 'regno')
        self.regnos = None
        self.editregno_btn.setEnabled(False)
        self.editregno_btn.clicked.connect(self.gotoEditRegno)
        self.editmol_btn.clicked.connect(self.editMolFile)
        self.regno_search_eb.editingFinished.connect(
            lambda: self.searchEvent(self.regno_search_eb.text(), 'regno'))

        self.loadsdf_btn.clicked.connect(self.gotoLoadSdf)
        self.back_btn.clicked.connect(self.previousRegno)
        self.forward_btn.clicked.connect(self.nextRegno)
        self.addmeta_btn.clicked.connect(self.gotoAddMeta)

        updateScreen(self)
        self.populated = True
        self.submitter_search_cb.currentTextChanged.connect(
            lambda x: self.searchEvent(x, 'CHEMIST'))

        self.compoundid_search_eb.textChanged.connect(
            lambda: self.searchEvent(self.compoundid_search_eb.text(), 'COMPOUND_ID'))

        self.batch_search_eb.editingFinished.connect(
            lambda: self.searchEvent(self.batch_search_eb.text(), 'JPAGE'))

    def editMolFile(self):
        fname = "tmp.mol" # temp file name
        fname_path = resource_path(fname) # file location / actual file name
        print(self.regno)
        if self.regno != None:
            # download molfile for selected regno, write to 'tmp.mol'
            tmp_mol_str = dbInterface.getMolFile(self.token,
                                                 self.regno,
                                                 self.database)
            tmp_file = open(fname_path, "w")
            n = tmp_file.write(tmp_mol_str)
            tmp_file.close()
            open_file(fname_path)

    def gotoLoadSdf(self):
        loadSDF = LoadSDF(self.token, self.database)

    def previousRegno(self):
        if self.regno == None or self.regnos == None:
            return
        self.searchingInProgress = True
        newIndex = self.regnos.index(self.regno) -1
        if newIndex > -1:
            self.regno = str(self.regnos[newIndex])
            self.numberOfHits_lab.setText(str(newIndex +1) + " of " + str(len(self.regnos)))
            self.populated == False
            updateScreen(self)
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
            updateScreen(self)
            self.populated == True
        self.searchingInProgress = False

    def gotoEditRegno(self):
        from regscreen import RegScreen
        reg = RegScreen(self.token, self.database, self.regno)
        self.window().addWidget(reg)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)
        
    def searchEvent(self, value='', action='doNothing'):
        if action == 'doNothing':
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
                                              self.token,
                                              self.database)
            for tmpRegno in tmpInts:
                self.regnos.append(str(tmpRegno))
            if len(self.regnos) > 0:
                sString = "1 of " + str(len(self.regnos))
                self.numberOfHits_lab.setText(sString)
                self.regno = str(self.regnos[0])
            else:
                sString = "0 of " + str(len(self.regnos))
                self.numberOfHits_lab.setText(sString)
                self.regno = None
            updateScreen(self)
            self.searchingInProgress = False

    def gotoReg(self):
        from regscreen import RegScreen
        reg = RegScreen(self.token, self.database)
        self.window().addWidget(reg)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

    def gotoAddMeta(self):
        addMetaTags = AddMetaTags(self.token)
