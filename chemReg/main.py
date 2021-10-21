import sys
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIntValidator, QImage, QPixmap
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QMainWindow, QFileDialog
import mysql
from mysql.connector import connect, Error
import requests
import json
import dbInterface

def updateScreen(self):
    #self.regno_eb.textChanged.connect(self.getMolfile)
    if self.populated == False:
        self.regno_eb.setText(self.regno)
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
    
    # Set current values for this regno
    if self.regno != None:
        currentBatch = dbInterface.getTextColumn(self.token, 'JPAGE', self.regno)
        self.batch_eb.setText(currentBatch)

        submitter = dbInterface.getTextColumn(self.token, 'chemist', self.regno)
        self.submitter_cb.setCurrentText(submitter)


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
        self.gotoReg(self.jwt_token)

    def gotoReg(self, token):
        reg = RegScreen(token)
        widget.addWidget(reg)
        widget.setCurrentIndex(widget.currentIndex() + 1)

    def gotoSearch(self, token):
        search = SearchScreen(token)
        widget.addWidget(search)
        widget.setCurrentIndex(widget.currentIndex() + 1)


class RegScreen(QMainWindow):
    def __init__(self, token, regnos=None):
        super(RegScreen, self).__init__()
        loadUi("regchem.ui", self)
        self.token = token
        self.dirty = False
        self.gotosearch_btn.clicked.connect(self.gotoSearch)
        self.setWindowTitle("Register new compound")
        self.onlyInt = QIntValidator()
        self.populated = False

        #### Regno
        if regnos == None:
            self.regno = dbInterface.getNextRegno(self.token)
            dbInterface.createNewRegno(self.regno, self.token)
        else:
            self.regno = regnos[0]
        self.regno_eb.textChanged.connect(self.getMolfile)
        self.regno_eb.setText(self.regno)

        updateScreen(self)
        self.populated = True
        #### Chemist
        #submitters = dbInterface.getColComboData(self.token,
        #                                         self.regno,
        #                                         'chemist')
        #self.submitter_cb.addItems(submitters)
        self.submitter_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'CHEMIST'))

        #### Projects
        #projects = dbInterface.getColComboData(self.token,
        #                                       self.regno,
        #                                       'project')
        #self.project_cb.addItems(projects)
        self.project_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'project'))

        #### CompoundTypes
        #compoundTypes = dbInterface.getColComboData(self.token,
        #                                            self.regno,
        #                                            'compound_type')
        #self.compoundtype_cb.addItems(compoundTypes)
        self.compoundtype_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'compound_type'))

        #### ProductTypes
        #productTypes = dbInterface.getColComboData(self.token,
        #                                           self.regno,
        #                                           'product')
        #self.product_cb.addItems(productTypes)
        self.product_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'product'))

        #### Libraries and LibraryName
        self.libraryid_cb.currentTextChanged.connect(self.changeLibraryName)
        #libraryIds = dbInterface.getColComboData(self.token,
        #                                         self.regno,
        #                                         'library_id')
        #self.libraryid_cb.addItems(libraryIds)
        self.libraryid_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'library_id'))

        #### ExternalId
        self.externalid_eb.editingFinished.connect(
            lambda: self.changeEvent(self.externalid_eb.text(), 'EXTERNAL_ID'))

        #### ExternalBatch
        self.externalbatch_eb.editingFinished.connect(
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
            
            sFile = "http://esox3.scilifelab.se:8082/mols/" + self.regno + ".png"
            image = QImage()
            self.structure_lab.setScaledContents(True)
            image.loadFromData(requests.get(sFile).content)
            self.structure_lab.setPixmap(QPixmap(image))
        
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
        self.regno = None
        self.regno_eb.textChanged.connect(
            lambda: self.searchEvent(self.regno_eb.text(), 'regno'))

        updateScreen(self)
        self.populated = True
        #submitters = dbInterface.getColComboData(self.token,
        #                                         None,
        #                                         'chemist')
        #self.submitter_cb.addItems(submitters)
        self.submitter_cb.currentTextChanged.connect(
            lambda x: self.searchEvent(x, 'CHEMIST'))

        self.compoundid_eb.textChanged.connect(
            lambda: self.searchEvent(self.compoundid_eb.text(), 'COMPOUND_ID'))

        self.batch_eb.textChanged.connect(
            lambda: self.searchEvent(self.batch_eb.text(), 'JPAGE'))


    def searchEvent(self, value='', action='doNothing'):
        if action == 'doNothing':
            return
        # Call search here
        regnos = dbInterface.searchValue(action, value, self.token)
        if len(regnos) > 0:
            self.regno = regnos[0]
            updateScreen(self)
            #reg = RegScreen(self.token, regnos)
            #widget.addWidget(reg)
            #widget.setCurrentIndex(widget.currentIndex() + 1)

    def gotoReg(self):
        reg = RegScreen(self.token)
        widget.addWidget(reg)
        widget.setCurrentIndex(widget.currentIndex() + 1)

        
app = QApplication(sys.argv)
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
