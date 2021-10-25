import sys
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIntValidator, QImage, QPixmap
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QMainWindow, QFileDialog
import mysql
from mysql.connector import connect, Error
import requests
import json
import dbInterface
import os

def displayMolfile(self):
    sFile = "http://esox3.scilifelab.se:8082/mols/" + self.regno + ".png"
    image = QImage()
    self.structure_lab.setScaledContents(True)
    image.loadFromData(requests.get(sFile).content)
    self.structure_lab.setPixmap(QPixmap(image))

def updateScreen(self):
    if self.populated == False:
        self.regno_eb.setText(self.regno)
        submitters = dbInterface.getColComboData(self.token, 'chemist')
        self.submitter_cb.addItems(submitters)
        
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
        
class LoginScreen(QDialog):
    def __init__(self):
        super(LoginScreen, self).__init__()
        loadUi("welcomescreen.ui", self)
        self.passwordfield.setEchoMode(QtWidgets.QLineEdit.Password)
        self.login.clicked.connect(self.loginfunction)
        
    def loginfunction(self):
        user = self.usernamefield.text()
        password = self.passwordfield.text()

        if len(user) == 0 or len(password) == 0:
            self.errorlabel.setText("Please input all fields")
        else:
            self.errorlabel.setText("")

        r = requests.post('http://esox3.scilifelab.se:8082/login',
                          data = {'username':user, 'password':password})
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
        loadUi("regchem.ui", self)
        self.token = token
        self.dirty = False
        self.gotosearch_btn.clicked.connect(self.gotoSearch)
        self.setWindowTitle("Register new compound")
        self.onlyInt = QIntValidator()
        self.populated = False

        #### Regno
        if regno == None:
            self.regno = dbInterface.getNextRegno(self.token)
            dbInterface.createNewRegno(self.regno, self.token)
        else:
            self.regno = regno
        self.regno_eb.textChanged.connect(self.getMolfile)
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
            f = {'file': open(fname[0], 'rb'),
                 'regno': self.regno}
            r = requests.post('http://esox3.scilifelab.se:8082/api/loadMolfile',
                              headers={'token': self.token}, files=f)
            displayMolfile(self)
        
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
    

class SearchScreen(QMainWindow):
    def __init__(self, token):
        super(SearchScreen, self).__init__()
        self.token = token
        loadUi("searchchem.ui", self)
        self.gotoreg_btn.clicked.connect(self.gotoReg)
        self.setWindowTitle("Search")
        self.populated = False
        self.searchingInProgress = False
        self.regno = None
        self.regnos = None
        self.editregno_btn.setEnabled(False)
        self.editregno_btn.clicked.connect(self.gotoEditRegno)
        self.regno_eb.editingFinished.connect(
            lambda: self.searchEvent(self.regno_eb.text(), 'regno'))

        self.back_btn.clicked.connect(self.previousRegno)
        self.forward_btn.clicked.connect(self.nextRegno)

        updateScreen(self)
        self.populated = True
        self.submitter_cb.currentTextChanged.connect(
            lambda x: self.searchEvent(x, 'CHEMIST'))

        self.compoundid_eb.textChanged.connect(
            lambda: self.searchEvent(self.compoundid_eb.text(), 'COMPOUND_ID'))

        self.batch_eb.editingFinished.connect(
            lambda: self.searchEvent(self.batch_eb.text(), 'JPAGE'))

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


# Handle high resolution displays:
if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"

app = QApplication(sys.argv)
app.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)

welcome = LoginScreen()
widget = QtWidgets.QStackedWidget()
widget.addWidget(welcome)
widget.setFixedHeight(800)
widget.setFixedWidth(1200)

widget.show()
try:
    sys.exit(app.exec_())
except:
    print("Exiting")
