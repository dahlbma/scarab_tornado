import sys, dbInterface, os, shutil, logging, re
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIntValidator
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from PyQt5.Qt import QApplication, QClipboard

from chemreglib import *


class EditStructure(QMainWindow):
    def __init__(self, token):
        super(EditStructure, self).__init__()
        loadUi(resource_path("assets/structure.ui"), self)
        self.token = token
        self.mod_name = "reg"
        logger = logging.getLogger(self.mod_name)
        self.window().setWindowTitle("Edit structure")
        self.compound_id = ''
        self.regno = ''
        
        self.compoundId_eb.textChanged.connect(self.compoundIdChanged)
        self.regno_eb.textChanged.connect(self.regnoChanged)

        self.editStructure_btn.clicked.connect(self.editMolFile)
        self.editStructure_btn.setEnabled(True)

        self.regnoCopyStructure_btn.clicked.connect(self.regnoCopyStructure)
        
        self.regnoForward_btn.clicked.connect(self.regnoForward)
        self.regnoBackward_btn.clicked.connect(self.regnoBackward)

        self.compoundForward_btn.clicked.connect(self.compoundForward)
        self.compoundBackward_btn.clicked.connect(self.compoundBackward)

        self.updateStructure_btn.setEnabled(False)
        self.updateStructure_btn.clicked.connect(self.updateStructure)


    def uploadMolfile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                            '.', "Molfiles (*.mol)")
        if fname[0] != '':
            #res, sMessage = postMolFile(self, fname[0], self.regno, logging.getLogger(self.mod_name))
            displayMolfile(self)
            updateMoleculeProperties(self)
            self.molChanged()      

    def editMolFile(self, event=None):
        self.fname = "tmp.mol" # temp file name
        self.fname_path = resource_path(self.fname) # file location / actual file name
        shutil.copy(resource_path("nostruct.mol"), self.fname_path)
        retcode = open_file(self.fname_path)

        self.fs_watcher = QtCore.QFileSystemWatcher([self.fname_path])
        self.fs_watcher.fileChanged.connect(self.createMolImage)        

    def createMolImage(self):
        molfile_file = open(self.fname_path, "r")
        molfile = molfile_file.read()
        molfile_file.close()
        self.statusMsg_lab.setText('')

        res = dbInterface.createMolImageFromMolfile(self.token, molfile)
        try:
            res = json.loads(res)
        except:
            print('Failed to create molecule image')
            return
            
        self.molfileToRegister = res['molfile']
        self.smilesToRegister = res['smiles']
        if res['status'] != '':
            self.statusMsg_lab.setText(res['status'])
        displayMolfile(self, 'new_structure', self.new_structure_lab)
        self.molfile = molfile
        self.updateStructure_btn.setEnabled(True)

    def compoundIdChanged(self):
        sCompoundId = self.compoundId_eb.text().upper()
        #pattern = '^CBK[0-9]{6}$'
        pattern = '^CBK[0-9]{6}'
        if len(re.findall(pattern, sCompoundId)) == 1:
            self.compoundId_eb.setText(sCompoundId)
            displayMolfile(self, 'no_struct', self.new_structure_lab)
            self.statusMsg_lab.setText('')
            self.compound_id = sCompoundId
            molkeyct = dbInterface.getMolkeyct(self.token, sCompoundId)
            self.molkeyct_eb.setText(str(molkeyct))
            saDuplicates = dbInterface.getCompoundDuplicates(self.token, sCompoundId)
            sConcatedCmps = ''
            for element in saDuplicates:
                sConcatedCmps += ' ' + element[0]

            self.cmpDuplicates_eb.setText(sConcatedCmps.strip())
            res = dbInterface.createBcpvsMolImage(self.token, sCompoundId)
            displayMolfile(self, sCompoundId, self.current_structure_lab)
            res = dbInterface.getRegnoFromCompound(self.token, sCompoundId)
            if res != '':
                self.regno_eb.setText(res)
            else:
                self.regno_eb.setText('')
                displayMolfile(self, 'no_struct', self.original_structure_lab)
        else:
            self.updateStructure_btn.setEnabled(False)


    def regnoCopyStructure(self):
        sRegno = self.regno_eb.text()
        sSql = f'''select molfile from chem_reg.chem_info
                   where regno = '{sRegno}' '''
        sMolfile = dbInterface.getMolFile(self.token, sRegno)
        
        cb = QApplication.clipboard()
        cb.clear(mode=cb.Clipboard)
        cb.setText(sMolfile, mode=cb.Clipboard)
        
        
    def compoundForward(self):
        sCmpId = self.compoundId_eb.text()
        forwardCompound = dbInterface.getForwardCompound(self.token, sCmpId)
        forwardCompound = forwardCompound.replace('"', '')
        if forwardCompound != '':
            self.compoundId_eb.setText(forwardCompound)
            sRegno = dbInterface.getRegnoFromCompound(self.token, forwardCompound)
            if sRegno != '':
                self.regno_eb.setText(sRegno)
            else:
                #self.compoundId_eb.setText('')
                displayMolfile(self, 'no_struct', self.original_structure_lab)

            
    def compoundBackward(self):
        sCmpId = self.compoundId_eb.text()
        backwardCompound = dbInterface.getBackwardCompound(self.token, sCmpId)
        backwardCompound = backwardCompound.replace('"', '')
        if backwardCompound != '':
            self.compoundId_eb.setText(backwardCompound)
            sRegno = dbInterface.getRegnoFromCompound(self.token, backwardCompound)
            if sRegno != '':
                self.regno_eb.setText(sRegno)
            else:
                #self.compoundId_eb.setText('')
                displayMolfile(self, 'no_struct', self.original_structure_lab)

            
    def regnoForward(self):
        sRegno = self.regno_eb.text()
        forwardRegno = dbInterface.getForwardRegno(self.token, sRegno)
        if forwardRegno != '':
            self.regno_eb.setText(forwardRegno)
            sCmpId = dbInterface.getCompoundFromRegno(self.token, forwardRegno)
            if sCmpId != '':
                self.compoundId_eb.setText(sCmpId)
            else:
                self.regno_eb.setText('')
                displayMolfile(self, 'no_struct', self.current_structure_lab)

        
    def regnoBackward(self):
        sRegno = self.regno_eb.text()
        prevRegno = dbInterface.getBackwardRegno(self.token, sRegno)
        if prevRegno != '':
            self.regno_eb.setText(prevRegno)
            sCmpId = dbInterface.getCompoundFromRegno(self.token, prevRegno)
            if sCmpId not in  ('', 'null'):
                self.compoundId_eb.setText(sCmpId)
            else:
                self.regno_eb.setText('')
                displayMolfile(self, 'no_struct', self.current_structure_lab)

        
    def regnoChanged(self):
        sRegno = self.regno_eb.text()
        #3913654
        #pattern = '^[0-9]{7}'
        pattern = '^[0-9]'
        if len(re.findall(pattern, sRegno)) == 1:
            displayMolfile(self, 'no_struct', self.new_structure_lab)
            self.statusMsg_lab.setText('')
            self.regno = sRegno
            res = dbInterface.createMolImage(self.token, sRegno)
            displayMolfile(self, sRegno, self.original_structure_lab)
            sCmpId = dbInterface.getCompoundFromRegno(self.token, sRegno)
            if sCmpId != '':
                self.compoundId_eb.setText(sCmpId)
            else:
                #self.regno_eb.setText('')
                displayMolfile(self, 'no_struct', self.original_structure_lab)
        else:
            self.updateStructure_btn.setEnabled(False)


    def molChanged(self):
        if (not self.new_structure_lab.pixmap().isNull()):
            self.molOK = True
        else:
            self.molOK = False
            self.molOK = True
        
    def updateStructure(self):
        # Update the structure in db here
        #dbInterface.updateStructure(self.token, self.compound_id, self.molfile)
        print('Update structure through dbInterface')
        prevRegno = dbInterface.updateStructureAdmin(self.token,
                                                     self.compound_id,
                                                     self.molfileToRegister,
                                                     self.smilesToRegister)
        res = dbInterface.createBcpvsMolImage(self.token, self.compound_id)
        displayMolfile(self,
                       self.compound_id,
                       self.current_structure_lab)


