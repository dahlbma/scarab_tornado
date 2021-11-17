import sys, requests, json, dbInterface, os, subprocess, platform, shutil, datetime, traceback, logging
from PyQt5.QtGui import QImage, QPixmap
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QMessageBox


def send_msg(title, text, icon=QMessageBox.Information, e=None):
    msg = QMessageBox()
    msg.setWindowTitle(title)
    msg.setIcon(icon)
    msg.setText(text)
    if e is not None:
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
    # opens folders too
    if platform.system() == 'Windows':
        if filename == ".":
            proc = os.startfile(filename)
        else:
            proc = subprocess.Popen("start " + filename, shell=True)
    elif platform.system() == 'Linux':
        proc = subprocess.Popen("xdg-open " + filename, shell=True)
    elif platform.system() == 'Darwin':
        proc = subprocess.Popen("open " + filename, shell=True)
    return proc

def displayMolfile(self):
    dbInterface.createMolImage(self.token, self.regno)
    sFile = "http://esox3.scilifelab.se:8082/mols/" + self.regno + ".png"
    image = QImage()
    self.structure_lab.setScaledContents(True)
    image.loadFromData(requests.get(sFile).content)
    self.structure_lab.setPixmap(QPixmap(image))

def postMolFile(self, fname, regno, logger):
    logger.info("posting file %s to server", fname)
    f = {'file': open(fname, 'rb'), 'regno': regno}
    r = requests.post('http://esox3.scilifelab.se:8082/api/loadMolfile',
                      headers={'token': self.token}, files=f)

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

        commentsText = dbInterface.getTextColumn(self.token, 'comments', self.regno)
        self.comments_text.setPlainText(commentsText)

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
      
