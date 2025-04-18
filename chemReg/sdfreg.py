import dbInterface, os, datetime, re, codecs, logging
from PyQt5.uic import loadUi
from PyQt5 import QtCore
from PyQt5.QtWidgets import QDialog, QApplication
from PyQt5.QtWidgets import QFileDialog, QProgressBar, QMessageBox
from chemreglib import *
import chardet
from bs4 import UnicodeDammit
import openpyxl

ip_rights_list = [None, 'External rights', 'LCBKI', 'Commercial']

class LoadSDF(QDialog):
    def __init__(self, token):
        super(LoadSDF, self).__init__()
        self.token = token
        self.mod_name = "sdfreg"
        logger = logging.getLogger(self.mod_name)
        self.sdfilename = None
        self.nostruct_name = None
        self.nostructs = None
        self.iMolCount = 0
        self.iNrElnIds = None
        self.saElnIds = None
        self.ElnIdsOK = False
        self.iFreeElnSpace = 0
        loadUi(resource_path("assets/sdfReg.ui"), self)
        self.pbar.setValue(0)
        self.pbar.hide()
        submitters = dbInterface.getColComboData(self.token,
                                                 'chemist')
        self.submitter_cb.addItems(submitters)
        self.submitter_cb.currentTextChanged.connect(self.check_fields)

        compoundTypes = dbInterface.getColComboData(self.token,
                                                    'compound_type')
        self.compoundtype_cb.addItems(compoundTypes)
        self.compoundtype_cb.currentTextChanged.connect(self.check_fields)

        projects = dbInterface.getColComboData(self.token,
                                               'project')
        self.project_cb.addItems(projects)
        self.project_cb.currentTextChanged.connect(self.check_fields)

        suppliers = dbInterface.getColComboData(self.token,
                                                'supplier')
        self.supplier_cb.addItems(suppliers)
        self.supplier_cb.currentTextChanged.connect(self.check_fields)

        solvents = dbInterface.getColComboData(self.token,
                                               'solvent')
        self.solvent_cb.addItems(solvents)
        self.solvent_cb.currentTextChanged.connect(self.check_fields)

        productTypes = dbInterface.getColComboData(self.token,
                                                   'product')
        self.producttype_cb.addItems(productTypes)
        self.producttype_cb.currentTextChanged.connect(self.check_fields)

        libraryIds = dbInterface.getColComboData(self.token,
                                                 'library_id')
        self.library_cb.addItems(libraryIds)
        self.library_cb.currentTextChanged.connect(self.update_librarydesc)
        self.library_cb.currentTextChanged.connect(self.check_fields)

        self.librarydesc_eb.setText(None)

        self.upload_btn.setEnabled(False)
        self.upload_btn.clicked.connect(self.uploadSDFile)
        self.selectsdf_btn.clicked.connect(self.getSDFile)
        self.nostruct_btn.clicked.connect(self.getNostructFile)
        self.cancel_btn.clicked.connect(self.closeWindow)

        self.elnids_text.textChanged.connect(self.parseElnIds)
        self.sdfname_lab.setText(' ')

        self.ip_rights_cb.addItems(ip_rights_list)
        self.ip_rights_cb.currentTextChanged.connect(self.check_fields)

        self.exec_()
        self.show()

    def check_fields(self):
        if self.sdfilename == None or \
           self.submitter_cb.currentText() == '' or \
           self.compoundtype_cb.currentText() == '' or \
           self.project_cb.currentText() == '' or \
           self.supplier_cb.currentText() == '' or \
           self.solvent_cb.currentText() == '' or \
           self.producttype_cb.currentText() == '' or \
           self.library_cb.currentText() in ('', ' ') or \
           self.ElnIdsOK == False or \
           self.ip_rights_cb.currentText() == '' or \
           self.iMolCount >= self.iFreeElnSpace:
            self.upload_btn.setEnabled(False)
        else:
            self.upload_btn.setEnabled(True)

        if self.nostruct_name != None:
            if self.nostruct_name == None or \
               self.submitter_cb.currentText() == '' or \
               self.compoundtype_cb.currentText() == '' or \
               self.project_cb.currentText() == '' or \
               self.supplier_cb.currentText() == '' or \
               self.solvent_cb.currentText() == '' or \
               self.producttype_cb.currentText() == '' or \
               self.library_cb.currentText() in ('', ' ') or \
               self.ElnIdsOK == False or \
               self.ip_rights_cb.currentText() == '' or \
               self.iMolCount >= self.iFreeElnSpace:
                self.upload_btn.setEnabled(False)            
            else:
                self.upload_btn.setEnabled(True)
                
        if self.nostruct_name != None and self.sdfilename != None:
            self.upload_btn.setEnabled(False)            
            

    
            

    def update_librarydesc(self):
        library_name = dbInterface.getLibraryName(self.token,
                                                  self.library_cb.currentText())
        self.librarydesc_eb.setText(library_name)

    def parseElnIds(self):
        sIds = self.elnids_text.toPlainText()
        saStrings = sIds.split(' ')
        iElnIdsFound = 0
        pattern = '^[a-zA-Z0-9]{6}$'
        iFreeElnPages = 0
        for sId in saStrings:
            if len(re.findall(pattern, sId)) == 1:
                iElnIdsFound += 1
                iLastPage = dbInterface.getLastBatchOfEln(self.token, sId)
                iFreeElnPages += 1000 - iLastPage
                self.iFreeElnSpace = iFreeElnPages
        saStrings = list(set(saStrings))
        if iElnIdsFound == self.iNrElnIds and len(saStrings) == iElnIdsFound:
            self.saElnIds = saStrings
            self.ElnIdsOK = True
        else:
            self.saElnIds = None
            self.ElnIdsOK = False
        self.check_fields()

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

            # RDKit can't handle empty first line in molfile
            if iCount == 1:
                line = b'id' + line

            sMol += line
            if b'$$$$' in line:
                sMol = sMol.decode(errors='replace')
                return sMol
        return ""

    def to_bytes(self, s):
        if type(s) is bytes:
            return s
        elif type(s) is str:
            return codecs.encode(s, 'utf-8')
        else:
            raise TypeError("Expected bytes or string, but got %s." % type(s))

    def getTags(self, sMol):
        pattern = b'>\s*<(.+)>.*\n(.*)\n'
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
                dammit = UnicodeDammit(i[1])
                dValues['external_id'] = dammit.unicode_markup
            elif i[0] == str.encode(self.batchfield_cb.currentText()):
                dValues['supplier_batch'] = i[1]
            elif i[0] == str.encode(self.purity_cb.currentText()):
                dValues['purity'] = i[1]
        return dValues

    def uploadNostructs(self):
        iElnId = 0
        sCurrentEln = self.saElnIds[iElnId]
        iBatchCount = dbInterface.getLastBatchOfEln(self.token, sCurrentEln)

        def getTags(external_id, supplier_batch, mw, restriction_comment):
            if restriction_comment == None:
                restriction_comment = ''
            if mw == None:
                mw = 0
            dTags = {
                "external_id": external_id,
                "supplier_batch": supplier_batch,
                "purity": -1,
                "mw": mw,
                "restriction_comment": restriction_comment
            }
            dTags['chemist'] = self.submitter_cb.currentText()
            dTags['compound_type'] = self.compoundtype_cb.currentText()
            dTags['project'] = self.project_cb.currentText()
            dTags['source'] = self.supplier_cb.currentText()
            dTags['solvent'] = self.solvent_cb.currentText()
            dTags['product'] = self.producttype_cb.currentText()
            dTags['library_id'] = self.library_cb.currentText()
            dTags['ip_rights'] = self.ip_rights_cb.currentText()
            dTags['sdfile_sequence'] = 0
            return dTags
            
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        for nostruct in self.nostructs:

            iBatchCount += 1
            if iBatchCount == 1000:
                iElnId += 1
                sCurrentEln = self.saElnIds[iElnId]
                iBatchCount = dbInterface.getLastBatchOfEln(self.token, sCurrentEln)
                if iBatchCount == 0:
                    iBatchCount = 1
    
            dTags = getTags(nostruct[0], nostruct[1], nostruct[2], nostruct[3])
            dTags['jpage'] = sCurrentEln + str(iBatchCount).zfill(3)
            dbInterface.addNostructMol(dTags, self.token)
        QApplication.restoreOverrideCursor()
        send_msg("All done", f"Done")

        
    def uploadSDFile(self):
        
        if self.nostruct_name != None and self.nostructs != None:
            self.uploadNostructs()
            return

        mol_info = {'external_id': self.cmpidfield_cb.currentText()}
        f = open(self.sdfilename, "rb")
        currtime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        f_err_path = f"error_{currtime}.sdf"
        f_err_msg_path = f"error_msg_{currtime}.log"
        f_err = open(f_err_path, "wb")
        f_err_msg = open(f_err_msg_path, "w")
        lError = False
        iTickCount = 0
        iTicks = int(self.iMolCount / 97)
        progress = 0
        iElnId = 0
        iNewMols = 0
        iErrorMols = 0
        self.pbar.show()
        self.pbar.setValue(progress)
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        sCurrentEln = self.saElnIds[iElnId]
        iBatchCount = dbInterface.getLastBatchOfEln(self.token,
                                                    sCurrentEln)
        iSdfSequence = dbInterface.getSdfSequence(self.token)
        self.event_lab.setText('Register in ChemReg')
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
            #if lTags == [] or sMol == "":
            if len(sMol) < 4:
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
                iElnId += 1
                sCurrentEln = self.saElnIds[iElnId]
                iBatchCount = dbInterface.getLastBatchOfEln(self.token,
                                                            sCurrentEln)
                if iBatchCount == 0:
                    iBatchCount = 1

            dTags['jpage'] = sCurrentEln + str(iBatchCount).zfill(3)
            dTags['molfile'] = sMol.decode('latin-1')
            dTags['chemist'] = self.submitter_cb.currentText()
            dTags['compound_type'] = self.compoundtype_cb.currentText()
            dTags['project'] = self.project_cb.currentText()
            dTags['source'] = self.supplier_cb.currentText()
            dTags['solvent'] = self.solvent_cb.currentText()
            dTags['product'] = self.producttype_cb.currentText()
            dTags['library_id'] = self.library_cb.currentText()
            dTags['ip_rights'] = self.ip_rights_cb.currentText()
            dTags['sdfile_sequence'] = iSdfSequence

            lStatus, sMessage = dbInterface.chemRegAddMolFile(dTags,
                                                              self.token)
            if sMessage == b'newMolecule':
                iNewMols += 1
            if lStatus != True:
                iErrorMols += 1
                f_err.write(sMol)
                f_err_msg.write(f"{str(dTags['external_id'])} {str(sMessage)}\n")
                f_err_msg.flush()
                lError = True
        saRegnos = dbInterface.getRegnosFromSdfSequence(self.token,
                                                        iSdfSequence)
        self.event_lab.setText('Register in Compound database')
        iTickCount = 0
        progress = 0
        for sReg in saRegnos:
            iTickCount += 1
            if iTickCount == iTicks:
                progress += 1
                iTickCount = 0
                self.pbar.setValue(progress)
            QApplication.processEvents()
            dbInterface.bcpvsRegCompound(self.token,
                                         sReg)
        self.event_lab.setText('')
        QApplication.processEvents()
        self.pbar.hide()
        QApplication.restoreOverrideCursor()

        if lError == False:
            send_msg("SDFile upload done", f'''Uploaded {self.iMolCount} compounds,
 {iNewMols} new compounds, {self.iMolCount - iNewMols} old compounds''')
            # remove error files
            f_err.close()
            f_err_msg.close()
            os.remove(f_err_path)
            os.remove(f_err_msg_path)
        else:
            send_msg("Some errors occured",
                     f'''Failed to register {iErrorMols} molecules,
 see {f_err_msg_path} and {f_err_path},
 {iNewMols} new compounds, {self.iMolCount - iNewMols} old compounds''',
                     QMessageBox.Warning)
            f_err.close()
            f_err_msg.close()
            logging.getLogger(self.mod_name).error(f'''Wrote failed molecule registration to
 {f_err_path} and info to {f_err_msg_path}''')
            open_file(".")

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
                self.iNrElnIds = int((iMolCount / 999) + 1)
                self.nrofelnids_lab.setText(str(self.iNrElnIds))
        f.close()
        f = open(fname[0], "rb")
        sMol = b''
        for line in f:
            sMol += line
            if b'$$$$' in line:
                break
        sMol = sMol.decode(errors='replace')
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
        if len(self.sdfilename) > 40:
            self.sdfname_lab.setText(self.sdfilename[:5] + "..." + self.sdfilename[-15:])
        else:
            self.sdfname_lab.setText(self.sdfilename)


    def getNostructFile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                                '.', "Excel file (*.xlsx)")
        if fname[0] == '':
            return

        #f = open(fname[0], "rb")
        iMolCount = 0
        self.nostruct_name = fname[0]
        if len(self.nostruct_name) > 40:
            self.nostructName_lab.setText(self.nostruct_name[:5] + "..." + self.nostruct_name[-19:])
        else:
            self.nostructName_lab.setText(self.nostruct_name)

        try:
            workbook = openpyxl.load_workbook(fname[0])
            sheet = workbook.active  # Get the active sheet
            header_row = [cell.value for cell in sheet[1]]  # Read the first row (header)
            expected_columns = ['External ID', 'External Batch', 'MW', 'restriction_comment']

            if not all(col in header_row for col in expected_columns) or len(expected_columns) != len(header_row):
                missing_cols = set(expected_columns) - set(header_row)
                extra_cols = set(header_row) - set(expected_columns)
                sError = ''
                if missing_cols:
                    sError = f'''Error: Missing columns: {missing_cols}\nName columns: 'External ID', 'External Batch', 'MW', 'restriction_comment'  '''
                if extra_cols:
                    sError = f'''Error: Unexpected columns: {extra_cols}\nName columns: 'External ID', 'External Batch', 'MW', 'restriction_comment' '''
                if len(expected_columns) != len(header_row):
                    sError = f'''Error: Incorrect number of columns. Expected {len(expected_columns)}, found {len(header_row)}\n\
Name columns: 'External ID', 'External Batch', 'MW', 'restriction_comment' '''

                send_msg("Format error", sError)
                return None  # Return None if the columns are not found

            data = []
            iMolCount = 0
            for row in sheet.iter_rows(min_row=2):  # Iterate from the second row (skipping the header)
                iMolCount += 1
                self.iMolCount = iMolCount
                self.compoundcount_lab.setText(str(iMolCount))
                self.iNrElnIds = int((iMolCount / 999) + 1)
                self.nrofelnids_lab.setText(str(self.iNrElnIds))
                row_values = [cell.value for cell in row]
                data.append(row_values)

            self.nostructs = data

        except FileNotFoundError:
            print(f"Error: File not found: {filepath}")
            return None
        except openpyxl.utils.exceptions.InvalidFileException:
            #print(f"Error: Invalid Excel file: {filepath}")
            return None
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return None

