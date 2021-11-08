import sys
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIntValidator, QImage, QPixmap
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QMainWindow
from PyQt5.QtWidgets import QFileDialog, QProgressBar, QMessageBox
import requests
import json
import dbInterface
import os
import subprocess, platform
import re
import codecs
import traceback
import logging

from PyQt5 import QtGui, uic

def setLogger(name='logger',
                level=logging.DEBUG,
                file=os.path.join(".","chemreg.log")):
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # console logging
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)

    # file logging
    fh = logging.FileHandler(file)
    fh.setLevel(level)
    formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(message)s', datefmt='%m/%d/%Y %H:%M:%S')
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger

def error_handler(etype, value, tb):
    err_msg = "".join(traceback.format_exception(etype, value, tb))
    #print("error caught")
    logger.exception(err_msg + " from error handler\n")

def send_msg(title, text, icon=QMessageBox.Information, e=None):
    msg = QMessageBox()
    msg.setWindowTitle(title)
    msg.setIcon(icon)
    msg.setText(text)
    if e is not None: #  icon is QMessageBox.Warning:
        # add clipboard btn
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Save)
        buttonS = msg.button(QMessageBox.Save)
        buttonS.setText('Save to clipboard')
    msg.exec_()
    if e is not None:
        if msg.clickedButton() == buttonS:
            # copy to clipboard if clipboard button was clicked
            clipboard.setText(text)
            cb_msg = QMessageBox()
            cb_msg.setText(clipboard.text()+" \n\ncopied to clipboard!")
            cb_msg.exec_()

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def open_file(filename):
    # open file with default OS application
    if platform.system() == 'Windows':
        proc = subprocess.Popen("start " + filename, shell=True)
    else:
        proc = subprocess.Popen("xdg-open " + filename, shell=True)
    return proc

def displayMolfile(self):
    sFile = "http://esox3.scilifelab.se:8082/mols/" + self.regno + ".png"
    image = QImage()
    self.structure_lab.setScaledContents(True)
    image.loadFromData(requests.get(sFile).content)
    self.structure_lab.setPixmap(QPixmap(image))

def postMolFile(self, fname, regno):
    logger.info("posting file %s to server", fname)
    f = {'file': open(fname, 'rb'), 'regno': regno}
    r = requests.post('http://esox3.scilifelab.se:8082/api/loadMolfile', headers={'token': self.token}, files=f)


def updateMoleculeProperties(self):
    avgMolMass = dbInterface.getTextColumn(self.token,
                                           'C_MW',
                                           self.regno)
    self.avgmolmass_lab.setText(avgMolMass)
    molFormula = dbInterface.getTextColumn(self.token,
                                           'C_MF',
                                           self.regno)
    self.mf_eb.setText(molFormula)
    sCHNS = dbInterface.getTextColumn(self.token,
                                      'C_CHNS',
                                      self.regno)
    self.chns_eb.setText(sCHNS)
    monoIsoMass = dbInterface.getTextColumn(self.token,
                                            'C_MONOISO',
                                            self.regno)
    self.monoisomass_lab.setText(monoIsoMass)
    
    
def updateScreen(self):
    if self.populated == False:
        self.regno_eb.setText(self.regno)
        submitters = dbInterface.getColComboData(self.token, 'chemist')
        self.submitter_cb.addItems(submitters)
        try:
            self.submitter_search_cb.addItems(submitters)
        except:
            pass
        projects = dbInterface.getColComboData(self.token, 'project')
        self.project_cb.addItems(projects)

        compoundTypes = dbInterface.getColComboData(self.token, 'compound_type')
        self.compoundtype_cb.addItems(compoundTypes)

        productTypes = dbInterface.getColComboData(self.token, 'product')
        self.product_cb.addItems(productTypes)
        
        libraryIds = dbInterface.getColComboData(self.token, 'library_id')
        self.libraryid_cb.addItems(libraryIds)
    
    # Set current values for this regno
    if self.regno != None:
        try:
            self.editregno_btn.setEnabled(True)
        except:
            pass
        self.regno_eb.setText(self.regno)
        
        createdDate = dbInterface.getTextColumn(self.token, 'rdate', self.regno)
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

        currentBatch = dbInterface.getTextColumn(self.token, 'JPAGE', self.regno)
        self.batch_eb.setText(currentBatch)

        submitter = dbInterface.getTextColumn(self.token, 'chemist', self.regno)
        self.submitter_cb.setCurrentText(submitter)
        try:
            self.submitter_search_cb.setCurrentText(submitter)
        except:
            pass

        project = dbInterface.getTextColumn(self.token, 'project', self.regno)
        self.project_cb.setCurrentText(project)

        product = dbInterface.getTextColumn(self.token, 'product', self.regno)
        self.product_cb.setCurrentText(product)

        compoundType = dbInterface.getTextColumn(self.token,
                                                 'compound_type',
                                                 self.regno)
        self.compoundtype_cb.setCurrentText(compoundType)

        libraryId = dbInterface.getTextColumn(self.token, 'library_id', self.regno)
        self.libraryid_cb.setCurrentText(libraryId)

        library_name = dbInterface.getLibraryName(self.token,
                                                  self.libraryid_cb.currentText())
        self.librarydesc_eb.setText(library_name)

        chromText = dbInterface.getTextColumn(self.token, 'chrom_text', self.regno)
        self.chrom_text.setPlainText(chromText)
        
        nmrText = dbInterface.getTextColumn(self.token, 'nmr_text', self.regno)
        self.nmr_text.setPlainText(nmrText)
        
        msText = dbInterface.getTextColumn(self.token, 'ms_text', self.regno)
        self.ms_text.setPlainText(msText)

        chromPurity = dbInterface.getTextColumn(self.token, 'chrom_purity', self.regno)
        self.chrompurity_eb.setText(chromPurity)
        nmrPurity = dbInterface.getTextColumn(self.token, 'nmr_purity', self.regno)
        self.nmrpurity_eb.setText(chromPurity)
        msPurity = dbInterface.getTextColumn(self.token, 'ms_purity', self.regno)
        self.mspurity_eb.setText(msPurity)

        displayMolfile(self)
    else:
        self.editregno_btn.setEnabled(False)
        self.regno_eb.setText(None)
        self.compoundid_eb.setText(None)
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
        self.externalid_eb.setText(None)
        self.externalbatch_eb.setText(None)
        self.iupac_eb.setText(None)
        self.monoisomass_lab.setText(None)
        self.avgmolmass_lab.setText(None)
        self.date_lab.setText(None)
        self.structure_lab.clear()
        # search buttons
        try:
            self.regno_search_eb.setText(None)
            self.compoundid_search_eb.setText(None)
            self.batch_search_eb.setText(None)
            self.submitter_search_cb.setCurrentText(' ')
        except:
            pass
        
class LoginScreen(QDialog):
    def __init__(self):
        super(LoginScreen, self).__init__()
        loadUi(resource_path("welcomescreen.ui"), self)
        self.passwordfield.setEchoMode(QtWidgets.QLineEdit.Password)
        self.login.clicked.connect(self.loginfunction)
        
    def loginfunction(self):
        user = self.usernamefield.text()
        password = self.passwordfield.text()

        if len(user) == 0 or len(password) == 0:
            self.errorlabel.setText("Please input all fields")
        else:
            self.errorlabel.setText("")
        
        try:
            r = requests.post('http://esox3.scilifelab.se:8082/login',
                              data = {'username':user, 'password':password})
        except Exception as e:
            self.errorlabel.setText("Bad Connection")
            send_msg("Error Message", str(e), QMessageBox.Warning, e)
            logger.exception(str(e))
            return
        if r.status_code != 200:
            self.errorlabel.setText("Wrong username/password")
            return
        self.jwt_token = r.content
        self.gotoSearch(self.jwt_token)

    def gotoReg(self, token):
        reg = RegScreen(token)
        widget.addWidget(reg)
        widget.setCurrentIndex(widget.currentIndex() + 1)

    def gotoSearch(self, token):
        search = SearchScreen(token)
        widget.addWidget(search)
        widget.setCurrentIndex(widget.currentIndex() + 1)


class RegScreen(QMainWindow):
    def __init__(self, token, regno=None):
        super(RegScreen, self).__init__()
        loadUi(resource_path("regchem.ui"), self)
        self.token = token
        self.dirty = False
        self.loadsdf_btn.clicked.connect(self.gotoLoadSdf)
        self.gotosearch_btn.clicked.connect(self.gotoSearch)
        widget.setWindowTitle("Register new compound")
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
        
        self.comments_text.textChanged.connect(
            lambda: self.changeEvent(self.comments_text.toPlainText(), 'COMMENTS'))

        self.loadmol_btn.clicked.connect(self.uploadMolfile)
        self.editmol_btn.clicked.connect(self.editMolFile)

    def changeLibraryName(self):
        library_name = dbInterface.getLibraryName(self.token,
                                                  self.libraryid_cb.currentText())
        self.librarydesc_eb.setText(library_name)

    def getMolfile(self):
        sFile = "http://esox3.scilifelab.se:8082/mols/" + self.regno + ".png"
        image = QImage()
        self.structure_lab.setScaledContents(True)
        image.loadFromData(requests.get(sFile).content)
        self.structure_lab.setPixmap(QPixmap(image))

    
        
    def uploadMolfile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                                '.', "Molfiles (*.mol)")
        if fname[0] != '':
            postMolFile(self, fname[0], self.regno)
            displayMolfile(self)
            updateMoleculeProperties(self)
            
    def editMolFile():
        fname = "sketch_test.mol"
        retcode = open_file(fname)
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
            #postMolFile(self, fname, self.regno)
            #displayMolfile(self)
            ok_msg.setText("Updated .mol file")
            print("y")
        else:
            # cancel, do nothing
            ok_msg.setText("Did not update .mol file")
        ok_msg.exec_()    
    
    def changeEvent(self, value='', action='doNothing'):
        if action == 'doNothing':
            return
        self.dirty = True
        dbInterface.updateValue(action, value, self.token, self.regno)

    def gotoSearch(self):
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno, self.token)
        search = SearchScreen(self.token)
        widget.addWidget(search)
        widget.setCurrentIndex(widget.currentIndex() + 1)

    def gotoLoadSdf(self):
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno, self.token)
        loadSDF = LoadSDF(self.token)


class LoadSDF(QDialog):
    def __init__(self, token):
        super(LoadSDF, self).__init__()
        self.token = token
        self.iMolCount = 0
        self.iNrElnIds = None
        self.saElnIds = None
        loadUi(resource_path("sdfReg.ui"), self)
        self.pbar.setValue(0)
        self.pbar.hide()
        submitters = dbInterface.getColComboData(self.token, 'chemist')
        self.submitter_cb.addItems(submitters)
        compoundTypes = dbInterface.getColComboData(self.token, 'compound_type')
        self.compoundtype_cb.addItems(compoundTypes)

        projects = dbInterface.getColComboData(self.token, 'project')
        self.project_cb.addItems(projects)

        suppliers = dbInterface.getColComboData(self.token, 'supplier')
        self.supplier_cb.addItems(suppliers)
        
        solvents = dbInterface.getColComboData(self.token, 'solvent')
        self.solvent_cb.addItems(solvents)
        
        productTypes = dbInterface.getColComboData(self.token, 'product')
        self.producttype_cb.addItems(productTypes)
        
        libraryIds = dbInterface.getColComboData(self.token, 'library_id')
        self.library_cb.addItems(libraryIds)

        self.upload_btn.setEnabled(False)
        self.upload_btn.clicked.connect(self.uploadSDFile)
        self.selectsdf_btn.clicked.connect(self.getSDFile)
        self.cancel_btn.clicked.connect(self.closeWindow)

        self.elnids_text.textChanged.connect(self.parseElnIds)

        self.exec_()
        self.show()

    def parseElnIds(self):
        sIds = self.elnids_text.toPlainText()
        saStrings = sIds.split(' ')
        iElnIdsFound = 0
        pattern = '^[a-zA-Z0-9]{9}$'
        for sId in saStrings:
            if len(re.findall(pattern, sId)) == 1:
                iElnIdsFound += 1
                if iElnIdsFound == self.iNrElnIds and \
                   len(saStrings) == self.iNrElnIds:
                    self.upload_btn.setEnabled(True)
                    self.saElnIds = saStrings
                else:    
                    self.upload_btn.setEnabled(False)

    def getNextMolecule(self, sFile):
        sMol = b""
        iCount = 0
        while True:
            iCount += 1
            if iCount > 1000:
                return ""
            line = sFile.readline()
            line = line.replace(b'\r\n', b'\n')
            line = line.replace(b"'", b"")
            sMol += line
            if b'$$$$' in line:
                return sMol[:-1] # + b"'"
        return ""

    def to_bytes(self, s):
        if type(s) is bytes:
            return s
        elif type(s) is str or type(s) is unicode:
            return codecs.encode(s, 'utf-8')
        else:
            raise TypeError("Expected bytes or string, but got %s." % type(s))

    
    def getTags(self, sMol):
        sPrevLine = ""
        pattern = b'>\s*<(.+)>\n(.*)\n'
        saTags = re.findall(pattern, sMol)
        return saTags
    
    def getValuePairs(self, lList):
        dValues = {
            "external_id": '',
            "supplier_batch": '',
            "purity": -1
            }

        for i in lList:            
            if i[0] == str.encode(self.cmpidfield_cb.currentText()):
                dValues['external_id'] = i[1]
            elif i[0] == str.encode(self.batchfield_cb.currentText()):
                dValues['supplier_batch'] = i[1]
            elif i[0] == str.encode(self.purity_cb.currentText()):
                dValues['purity'] = i[1]
        return dValues
        
    def uploadSDFile(self):
        mol_info = {'external_id': self.cmpidfield_cb.currentText()}
        f = open(self.sdfilename, "rb")
        f_err = open("error.sdf", "wb")
        f_err_msg = open("error_msg.txt", "w")
        lError = False
        iTickCount = 0
        iBatchCount = 0
        iTicks = int(self.iMolCount / 100)
        progress = 0
        iElnId = 0
        self.pbar.show()
        self.pbar.setValue(progress)
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        sCurrentEln = self.saElnIds[iElnId]
        while True:
            # Maybe run in separate thread instead, see:
            # https://www.pythonguis.com/tutorials/multithreading-pyqt-applications-qthreadpool/
            QApplication.processEvents()
            iTickCount += 1
            if iTickCount == iTicks:
                progress += 1
                iTickCount = 0
                self.pbar.setValue(progress)
            sMol = self.getNextMolecule(f)
            sMol = self.to_bytes(sMol)
            lTags = self.getTags(sMol)
            if lTags == [] or sMol == "":
                break
            dTags = self.getValuePairs(lTags)

            try:
                iSlask = int(dTags['purity']) + 1
            except:
                send_msg("Purity error", f"Purity must be a number, got: {dTags['purity']}")
                lError = True
                break
            
            iBatchCount += 1
            if iBatchCount == 1000:
                iBatchCount = 1
                iElnId += 1
                sCurrentEln = self.saElnIds[iElnId]
            dTags['jpage'] = sCurrentEln + str(iBatchCount).zfill(3)
            dTags['molfile'] = sMol.decode('latin-1')
            dTags['chemist'] = self.submitter_cb.currentText()
            dTags['compound_type'] = self.compoundtype_cb.currentText()
            dTags['project'] = self.project_cb.currentText()
            dTags['source'] = self.supplier_cb.currentText()
            dTags['solvent'] = self.solvent_cb.currentText()
            dTags['product'] = self.producttype_cb.currentText()
            dTags['library_id'] = self.library_cb.currentText()
            lStatus, sMessage = dbInterface.uploadMolFile(dTags, self.token)
            if lStatus != True:
                f_err.write(sMol)
                f_err_msg.write(f"{str(dTags['external_id'])} {str(sMessage)}\n")
                f_err_msg.flush()
                lError = True
        self.pbar.hide()
        QApplication.restoreOverrideCursor()
        if lError == False:
            send_msg("SDFile upload done", f"Uploaded {self.iMolCount} compounds")
        else:
            send_msg("Some errors occured",
                     f"Failed to register some molecules, see error_msg.txt and error.sdf",
                     QMessageBox.Warning)

            
    def closeWindow(self):
        self.close()
    
    def getSDFile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                                '.', "SDFiles (*.sdf)")
        if fname[0] == '':
            return

        f = open(fname[0], "rb")
        iMolCount = 0
        for line in f:
            if b'$$$$' in line:
                iMolCount += 1
                self.iMolCount = iMolCount
                self.compoundcount_lab.setText(str(iMolCount))
                self.iNrElnIds = int((iMolCount / 1000) + 1)
                self.nrofelnids_lab.setText(str(self.iNrElnIds))
        f.close()
        f = open(fname[0], "rb")
        sMol = b''
        for line in f:
            sMol += line
            if b'$$$$' in line:
                break
        sMol = sMol.decode()
        pattern = '>\s*<(.+)>'
        saTags = re.findall(pattern, sMol)
        saTags.insert(0, '')
        self.cmpidfield_cb.clear()
        self.cmpidfield_cb.addItems(saTags)
        self.batchfield_cb.clear()
        self.batchfield_cb.addItems(saTags)
        self.purity_cb.clear()
        self.purity_cb.addItems(saTags)
        self.sdfilename = fname[0]


class SearchScreen(QMainWindow):
    def __init__(self, token):
        super(SearchScreen, self).__init__()
        self.token = token
        loadUi(resource_path("searchchem.ui"), self)
        self.gotoreg_btn.clicked.connect(self.gotoReg)
        widget.setWindowTitle("Search")
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
        reg = RegScreen(self.token, self.regno)
        widget.addWidget(reg)
        widget.setCurrentIndex(widget.currentIndex() + 1)
        
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
        reg = RegScreen(self.token)
        widget.addWidget(reg)
        widget.setCurrentIndex(widget.currentIndex() + 1)


logger = setLogger()
sys.excepthook = error_handler

try:
    os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "2"
    app = QApplication(['Chem Reg'])
    clipboard = app.clipboard()

    welcome = LoginScreen()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(welcome)

    desktop = QApplication.desktop()
    windowHeight = 800
    windowWidth = 1200

    windowHeight = int(round(0.9 * desktop.screenGeometry().height(), -1))
    if windowHeight > 800:
        windowHeight = 800

    windowWidth = int(round((1200/800) * windowHeight, -1))

    widget.resize(windowWidth, windowHeight)

    widget.show()
    sys.exit(app.exec_())
except Exception as e:
    logger.info(str(e))
