import dbInterface, os, datetime, re, codecs, logging, math
from PyQt5.uic import loadUi
from PyQt5 import QtCore
from PyQt5.QtWidgets import QDialog, QApplication
from PyQt5.QtWidgets import QFileDialog, QProgressBar, QMessageBox
from chemreglib import *
import chardet
from bs4 import UnicodeDammit
import openpyxl
import pandas as pd

ip_rights_list = [None, 'External rights', 'LCBKI', 'Commercial']
NOSTRUCT_COLUMNS = ('External ID', 'External Batch', 'MW', 'restriction_comment')

class LoadSDF(QDialog):
    def __init__(self, token):
        super(LoadSDF, self).__init__()
        self.token = token
        self.mod_name = "sdfreg"
        logger = logging.getLogger(self.mod_name)
        self.sdfilename = None
        self.excelfilename = None
        self.df = None
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
        self.selectSmiles_btn.clicked.connect(self.getExcelSmileFile)
        self.nostruct_btn.clicked.connect(self.getNostructFile)
        self.cancel_btn.clicked.connect(self.closeWindow)

        self.elnids_text.textChanged.connect(self.parseElnIds)
        self.sdfname_lab.setText(' ')

        self.ip_rights_cb.addItems(ip_rights_list)
        self.ip_rights_cb.currentTextChanged.connect(self.check_fields)

        self.exec_()
        self.show()

    def check_fields(self):
        # Special handling for CTRL Compound - skip the iMolCount check
        is_ctrl_compound = self.compoundtype_cb.currentText() == 'CTRL Compound'
        
        if is_ctrl_compound:
            # Block signals to prevent recursion, then set the values
            self.elnids_text.blockSignals(True)
            self.elnids_text.setPlainText('CTRL')
            self.elnids_text.blockSignals(False)
            # Set the flags directly for CTRL compounds
            self.saElnIds = ['CTRL']
            self.ElnIdsOK = True
            self.iFreeElnSpace = 1000000
        else:
            # If switching away from CTRL Compound, clear the CTRL text and reset flags
            if self.elnids_text.toPlainText() == 'CTRL':
                self.elnids_text.blockSignals(True)
                self.elnids_text.setPlainText('')
                self.elnids_text.blockSignals(False)
                self.saElnIds = None
                self.ElnIdsOK = False
                self.iFreeElnSpace = 0

        has_sdf = self.sdfilename is not None
        has_smiles = self.excelfilename is not None and self.df is not None
        has_nostruct = self.nostruct_name is not None and self.nostructs is not None

        if sum((has_sdf, has_smiles, has_nostruct)) != 1:
            self.upload_btn.setEnabled(False)
            return

        fields_are_valid = self.submitter_cb.currentText() != '' and \
                           self.compoundtype_cb.currentText() != '' and \
                           self.project_cb.currentText() != '' and \
                           self.supplier_cb.currentText() != '' and \
                           self.solvent_cb.currentText() != '' and \
                           self.producttype_cb.currentText() != '' and \
                           self.library_cb.currentText() not in ('', ' ') and \
                           self.ElnIdsOK and \
                           self.ip_rights_cb.currentText() != ''
        enough_eln_space = is_ctrl_compound or self.iMolCount < self.iFreeElnSpace
        self.upload_btn.setEnabled(fields_are_valid and enough_eln_space)



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

        # Special case for CTRL, which means that no ELN IDs are used
        if sIds.strip() == 'CTRL' and self.compoundtype_cb.currentText() == 'CTRL Compound':
            self.saElnIds = ['CTRL']
            self.ElnIdsOK = True
            self.iFreeElnSpace = 1000000
            self.check_fields()
            return
        
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
            # NOTE: do not strip apostrophes from supplier data here.
            # SQL injection is now prevented server-side via
            # parameterized queries, and stripping apostrophes silently
            # corrupted supplier identifiers (e.g. names containing "'").

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

    def getMolfileTags(self, sMol):
        pattern = rb'>\s*<(.+)>.*\n(.*)\n'
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

    def getExcelValuePairs(self, row):
        dValues = {
            "external_id": '',
            "supplier_batch": '',
            "purity": -1
        }

        for field_name, value_key in (
                (self.cmpidfield_cb.currentText(), 'external_id'),
                (self.batchfield_cb.currentText(), 'supplier_batch'),
                (self.purity_cb.currentText(), 'purity')):
            if not field_name:
                continue
            value = row[field_name]
            if pd.isna(value):
                continue
            if value_key == 'purity':
                try:
                    numeric_value = float(value)
                    if numeric_value.is_integer():
                        value = int(numeric_value)
                except (TypeError, ValueError):
                    pass
            dValues[value_key] = str(value)
        return dValues


    def getTags(self, external_id, supplier_batch, mw, restriction_comment):
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

    def validateNostructData(self, rows):
        errors = []
        validated_rows = []

        if not rows:
            return ["The workbook must contain at least one data row."], []

        for row_number, row in enumerate(rows, start=2):
            if len(row) != len(NOSTRUCT_COLUMNS):
                errors.append(
                    f"Row {row_number}: expected {len(NOSTRUCT_COLUMNS)} "
                    f"values, found {len(row)}."
                )
                continue

            external_id, supplier_batch, mw, restriction_comment = row
            row_errors = []
            external_id = '' if external_id is None else str(external_id).strip()
            supplier_batch = '' if supplier_batch is None else str(supplier_batch).strip()

            if not external_id:
                row_errors.append(f"Row {row_number}: External ID is required.")
            elif len(external_id) > 200:
                row_errors.append(
                    f"Row {row_number}: External ID must be at most 200 characters."
                )

            if not supplier_batch:
                row_errors.append(f"Row {row_number}: External Batch is required.")
            elif len(supplier_batch) > 100:
                row_errors.append(
                    f"Row {row_number}: External Batch must be at most 100 characters."
                )

            if mw is None or (isinstance(mw, str) and not mw.strip()):
                mw = 0.0
            else:
                try:
                    mw = float(mw)
                except (TypeError, ValueError):
                    row_errors.append(f"Row {row_number}: MW must be a number.")
                else:
                    if not math.isfinite(mw):
                        row_errors.append(
                            f"Row {row_number}: MW must be a finite number."
                        )

            restriction_comment = '' if restriction_comment is None else str(
                restriction_comment
            ).strip()

            if row_errors:
                errors.extend(row_errors)
                continue

            validated_rows.append([
                external_id,
                supplier_batch,
                mw,
                restriction_comment
            ])

        return errors, validated_rows

    def showNostructValidationError(self, errors):
        max_errors = 20
        error_text = '\n'.join(errors[:max_errors])
        if len(errors) > max_errors:
            error_text += f"\n... and {len(errors) - max_errors} more errors."
        send_msg(
            "Nostruct file format error",
            "Expected columns in this exact order:\n"
            + ", ".join(NOSTRUCT_COLUMNS)
            + "\n\n"
            "Row requirements:\n"
            "- External ID: non-empty text, at most 200 characters\n"
            "- External Batch: non-empty text, at most 100 characters\n"
            "- MW: a finite number; blank is treated as 0\n"
            "- restriction_comment: text; blank is allowed\n\n"
            "Problems found:\n"
            + error_text,
            QMessageBox.Warning
        )


    def uploadControls(self):
        print('Uploading controls')
        sCurrentEln = self.saElnIds[0]
        iBatchCount = dbInterface.getLastCtrlComp(self.token, sCurrentEln)
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        sCtrls = ""
        for nostruct in self.nostructs:

            iBatchCount += 1
            dTags = self.getTags(nostruct[0], nostruct[1], nostruct[2], nostruct[3])
            dTags['jpage'] = sCurrentEln + str(iBatchCount).zfill(5)
            dTags['iCtrlComp'] = iBatchCount
            sCtrls += dTags['jpage'] + " "
            dbInterface.addCtrlMol(dTags, self.token)
        QApplication.restoreOverrideCursor()
        send_msg("All done", f"Done, added control compounds: {sCtrls}")


    def uploadNostructs(self):
        validation_errors, validated_rows = self.validateNostructData(self.nostructs)
        if validation_errors:
            self.showNostructValidationError(validation_errors)
            return
        self.nostructs = validated_rows

        iElnId = 0
        sCurrentEln = self.saElnIds[iElnId]

        if sCurrentEln == 'CTRL' and self.compoundtype_cb.currentText() == 'CTRL Compound':
            self.uploadControls()
            return

        iBatchCount = dbInterface.getLastBatchOfEln(self.token, sCurrentEln)
            
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        for nostruct in self.nostructs:

            iBatchCount += 1
            if iBatchCount == 1000:
                iElnId += 1
                sCurrentEln = self.saElnIds[iElnId]
                iBatchCount = dbInterface.getLastBatchOfEln(self.token, sCurrentEln)
                if iBatchCount == 0:
                    iBatchCount = 1

            dTags = self.getTags(nostruct[0], nostruct[1], nostruct[2], nostruct[3])
            dTags['jpage'] = sCurrentEln + str(iBatchCount).zfill(3)
            dbInterface.addNostructMol(dTags, self.token)
        QApplication.restoreOverrideCursor()
        send_msg("All done", f"Done")

        
    def uploadSDFile(self):
        if self.nostruct_name != None and self.nostructs != None:
            self.uploadNostructs()
            return

        is_smiles_upload = self.excelfilename is not None and self.df is not None
        if is_smiles_upload:
            excel_rows = iter(self.df.to_dict(orient='records'))
            f = None
        else:
            excel_rows = None
            f = open(self.sdfilename, "rb")
        currtime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        error_extension = 'smi' if is_smiles_upload else 'sdf'
        f_err_path = f"error_{currtime}.{error_extension}"
        f_err_msg_path = f"error_msg_{currtime}.log"
        f_err = open(f_err_path, "wb")
        f_err_msg = open(f_err_msg_path, "w")
        lError = False
        iTickCount = 0
        iTicks = max(1, int(self.iMolCount / 97))
        progress = 0
        iElnId = 0
        iNewMols = 0
        iErrorMols = 0
        iPostedMols = 0
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
            if is_smiles_upload:
                try:
                    row = next(excel_rows)
                except StopIteration:
                    break
                raw_structure = row.get('Smile')
                sStructure = '' if pd.isna(raw_structure) else str(raw_structure).strip()
                dTags = self.getExcelValuePairs(row)
                if not sStructure or sStructure.lower() == 'nan':
                    iErrorMols += 1
                    f_err_msg.write(f"{dTags['external_id']} Empty Smile\n")
                    f_err_msg.flush()
                    lError = True
                    continue
            else:
                sMol = self.getNextMolecule(f)
                sMol = self.to_bytes(sMol)
                lTags = self.getMolfileTags(sMol)
                if len(sMol) < 4:
                    break
                dTags = self.getValuePairs(lTags)
                sStructure = sMol

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
            if is_smiles_upload:
                dTags['smiles'] = sStructure
            else:
                dTags['molfile'] = sStructure.decode('latin-1')
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
            if lStatus is True:
                iPostedMols += 1
            # Server may report a tautomer warning while still
            # registering the row successfully (status 200, body of
            # the form b"newMolecule;warning=..."). Capture these so
            # the user can review without aborting the upload.
            try:
                sMessageText = sMessage.decode('utf-8', errors='replace') if isinstance(sMessage, bytes) else str(sMessage)
            except Exception:
                sMessageText = ''
            if sMessageText.split(';warning=', 1)[0] == 'newMolecule':
                iNewMols += 1
            if lStatus is True and ';warning=' in sMessageText:
                f_err_msg.write(f"{str(dTags['external_id'])} WARNING {sMessageText}\n")
                f_err_msg.flush()
                if sMessageText.split(';warning=', 1)[0] == 'newMolecule':
                    pass  # already counted above
            if lStatus != True:
                iErrorMols += 1
                if is_smiles_upload:
                    f_err.write(sStructure.encode('utf-8') + b'\n')
                else:
                    f_err.write(sStructure)
                f_err_msg.write(f"{str(dTags['external_id'])} {str(sMessage)}\n")
                f_err_msg.flush()
                lError = True
        saRegnos = dbInterface.getRegnosFromSdfSequence(self.token,
                                                        iSdfSequence)
        if len(saRegnos) != iPostedMols:
            iErrorMols += 1
            f_err_msg.write(
                f"Expected {iPostedMols} registered rows for sequence "
                f"{iSdfSequence}, found {len(saRegnos)}\n")
            f_err_msg.flush()
            lError = True
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
            bcpvs_status, bcpvs_message = dbInterface.bcpvsRegCompound(
                self.token, sReg, return_status=True)
            if not bcpvs_status:
                iErrorMols += 1
                f_err_msg.write(f"{sReg} {str(bcpvs_message)}\n")
                f_err_msg.flush()
                lError = True
        self.event_lab.setText('')
        QApplication.processEvents()
        self.pbar.hide()
        QApplication.restoreOverrideCursor()

        if lError == False:
            upload_type = 'SMILES' if is_smiles_upload else 'SDFile'
            send_msg(f"{upload_type} upload done", f'''Uploaded {self.iMolCount} compounds,
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

        self.excelfilename = None
        self.df = None
        self.nostruct_name = None
        self.nostructs = None
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
        pattern = r'>\s*<(.+)>'
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

    def getExcelSmileFile(self):
        fname = QFileDialog.getOpenFileName(
            self, "Open file", ".", "Excel Files (*.xlsx *.xls)"
        )
        if not fname[0]:
            return

        df = pd.read_excel(fname[0])
        if 'Smile' not in df.columns:
            send_msg("Format error", 'The Excel file must contain a column named "Smile".')
            return

        self.sdfilename = None
        self.nostruct_name = None
        self.nostructs = None
        df.columns = [str(column) for column in df.columns]

        # Calculate row count and ELN batch splits
        iMolCount = len(df)
        self.iMolCount = iMolCount
        self.compoundcount_lab.setText(str(iMolCount))
        self.iNrElnIds = int((iMolCount / 999) + 1)
        self.nrofelnids_lab.setText(str(self.iNrElnIds))

        # Extract dynamic column names and prepend empty mapping option
        columns = list(df.columns)
        columns.insert(0, "")

        # Populate mapping dropdowns
        self.cmpidfield_cb.clear()
        self.cmpidfield_cb.addItems(columns)
        self.batchfield_cb.clear()
        self.batchfield_cb.addItems(columns)
        self.purity_cb.clear()
        self.purity_cb.addItems(columns)

        # Store file path and update label string truncation
        self.excelfilename = fname[0]
        if len(self.excelfilename) > 40:
            self.smilesname_lab.setText(
                self.excelfilename[:5] + "..." + self.excelfilename[-15:]
            )
        else:
            self.smilesname_lab.setText(self.excelfilename)

        # Store DataFrame reference for downstream downstream processing/parsing
        self.df = df
        self.check_fields()

    def getNostructFile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                                '.', "Excel file (*.xlsx)")
        if fname[0] == '':
            return

        self.sdfilename = None
        self.excelfilename = None
        self.df = None
        self.nostructs = None
        self.nostruct_name = None
        self.nostructName_lab.setText(' ')
        self.iMolCount = 0
        self.iNrElnIds = None
        self.compoundcount_lab.setText('0')
        self.nrofelnids_lab.setText('0')
        self.check_fields()

        iMolCount = 0
        try:
            workbook = openpyxl.load_workbook(fname[0])
            sheet = workbook.active  # Get the active sheet
            header_row = [cell.value for cell in sheet[1]]  # Read the first row (header)
            if header_row != list(NOSTRUCT_COLUMNS):
                self.showNostructValidationError([
                    "Header row does not match the required columns in the "
                    "required order.",
                    f"Found: {header_row}"
                ])
                return None  # Return None if the columns are not found

            data = [
                [cell.value for cell in row]
                for row in sheet.iter_rows(min_row=2)
            ]
            validation_errors, validated_rows = self.validateNostructData(data)
            if validation_errors:
                self.showNostructValidationError(validation_errors)
                return None

            self.nostructs = validated_rows
            self.nostruct_name = fname[0]
            if len(self.nostruct_name) > 40:
                self.nostructName_lab.setText(self.nostruct_name[:5] + "..." + self.nostruct_name[-19:])
            else:
                self.nostructName_lab.setText(self.nostruct_name)
            self.iMolCount = len(self.nostructs)
            self.compoundcount_lab.setText(str(self.iMolCount))
            self.iNrElnIds = int((self.iMolCount / 999) + 1)
            self.nrofelnids_lab.setText(str(self.iNrElnIds))
            # Trigger check_fields to re-evaluate button state with the loaded file
            self.check_fields()

        except FileNotFoundError:
            send_msg(
                "File error",
                f"Could not find the selected file: {fname[0]}",
                QMessageBox.Warning
            )
            return None
        except openpyxl.utils.exceptions.InvalidFileException:
            send_msg(
                "Format error",
                "The selected file is not a valid Excel workbook.",
                QMessageBox.Warning
            )
            return None
        except Exception as e:
            send_msg(
                "File error",
                f"Could not read the selected workbook:\n{e}",
                QMessageBox.Warning
            )
            return None

