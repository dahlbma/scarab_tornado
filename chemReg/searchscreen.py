import sys, dbInterface, os
from PyQt5.uic import loadUi
from PyQt5.QtWidgets import QMainWindow

from chemreglib import *

from sdfreg import LoadSDF


class SearchScreen(QMainWindow):
    def __init__(self, token):
        super(SearchScreen, self).__init__()
        self.token = token
        self.logger = setLogger()
        loadUi(resource_path("searchchem.ui"), self)
        self.gotoreg_btn.clicked.connect(self.gotoReg)
        self.window().setWindowTitle("Search")
        self.dirty = False
        self.populated = False
        self.searchingInProgress = False
        self.regno = None
        self.regnos = None
        self.editregno_btn.setEnabled(False)
        self.editregno_btn.clicked.connect(self.gotoEditRegno)
        self.regno_search_eb.editingFinished.connect(
            lambda: self.searchEvent(self.regno_search_eb.text(), 'regno'))

        self.loadsdf_btn.clicked.connect(self.gotoLoadSdf)
        self.back_btn.clicked.connect(self.previousRegno)
        self.forward_btn.clicked.connect(self.nextRegno)

        updateScreen(self)
        self.populated = True
        self.submitter_search_cb.currentTextChanged.connect(
            lambda x: self.searchEvent(x, 'CHEMIST'))

        self.compoundid_search_eb.textChanged.connect(
            lambda: self.searchEvent(self.compoundid_search_eb.text(), 'COMPOUND_ID'))

        self.batch_search_eb.editingFinished.connect(
            lambda: self.searchEvent(self.batch_search_eb.text(), 'JPAGE'))

    def gotoLoadSdf(self):
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno, self.token)
        loadSDF = LoadSDF(self.token)

    def previousRegno(self):
        if self.regno == None:
            return
        self.searchingInProgress = True
        newIndex = self.regnos.index(self.regno) -1
        if newIndex > -1:
            self.regno = self.regnos[newIndex]
            self.numberOfHits_lab.setText(str(newIndex +1) + " of " + str(len(self.regnos)))
            updateScreen(self)
        self.searchingInProgress = False
        
    def nextRegno(self):
        if self.regno == None:
            return
        self.searchingInProgress = True
        newIndex = self.regnos.index(self.regno) +1
        if newIndex < len(self.regnos):
            self.numberOfHits_lab.setText(str(newIndex +1) + " of " + str(len(self.regnos)))
            self.regno = self.regnos[newIndex]
            updateScreen(self)
        self.searchingInProgress = False

    def gotoEditRegno(self):
        from regscreen import RegScreen
        reg = RegScreen(self.token, self.regno)
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
            self.regnos = dbInterface.searchValue(action, value, self.token)
            if len(self.regnos) > 0:
                sString = "1 of " + str(len(self.regnos))
                self.numberOfHits_lab.setText(sString)
                self.regno = self.regnos[0]
            else:
                sString = "0 of " + str(len(self.regnos))
                self.numberOfHits_lab.setText(sString)
                self.regno = None
            updateScreen(self)
            self.searchingInProgress = False

    def gotoReg(self):
        from regscreen import RegScreen
        reg = RegScreen(self.token)
        self.window().addWidget(reg)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)
