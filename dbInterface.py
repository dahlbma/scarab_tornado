import tornado.web
import json
import MySQLdb
from auth import jwtauth
import config
import logging
import rdkit.Chem
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
from molmass import Formula
from io import StringIO
import sys
import codecs
import re
import os
import mydb
import chembl_export
import zipfile
import random

logger = logging.getLogger(__name__)

db = mydb.disconnectSafeConnect()
cur = db.cursor()


def getDatabase(parent):
    data = parent.request.headers['Token']
    jsonData = json.loads(data)
    database = jsonData['database']
    if database == 'Live':
        return 'chem_reg', 'bcpvs'
    else:
        return 'chem_reg_test', 'bcpvs_test'

def res2json():
    result = [list(i) for i in cur.fetchall()]
    return json.dumps(result)


def checkUniqueStructure(inchiKey, bcpvsDB):
    """Look up an existing compound by InChIKey.

    Historically this matched on `smiles_std`, but a canonical SMILES is
    sensitive to the (non-deterministic) RDKit tautomer canonicalizer
    and to atom ordering, which caused the same compound to be
    registered twice across different runs. The `inchi_key` column on
    `bcpvs.compound` is indexed and is what should be used for identity.
    """
    if not inchiKey:
        return ()
    sSql = f"SELECT compound_id FROM {bcpvsDB}.`compound` WHERE inchi_key = %s"
    cur.execute(sSql, (inchiKey,))
    return cur.fetchall()


def getNewRegno(chemregDB):
    sSql = f"select pkey from {chemregDB}.regno_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update {chemregDB}.regno_sequence set pkey={pkey}"
    cur.execute(sSql)
    return pkey

def getNewSaltNumber(chemregDB):
    sSql = f"select pkey from {chemregDB}.salts_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update {chemregDB}.salts_sequence set pkey={pkey}"
    cur.execute(sSql)
    return pkey

def getNewCompoundId(bcpvsDB):
    sSql = f"select pkey from {bcpvsDB}.compound_id_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update {bcpvsDB}.compound_id_sequence set pkey={pkey}"
    cur.execute(sSql)
    return pkey

def getAtomicComposition(saComp):
    sComp = f""
    for atom in saComp:
        sComp += f'{atom[0]} {round(atom[3] * 100, 2)}% '
    return sComp

def createPngFromMolfile(regno, molfile):
    m = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    try:
        Draw.MolToFile(m, f'mols/{regno}.png', kekulize=True, size=(280, 280))
    except:
        logger.error(f"regno {regno} is nostruct")

def updateCtrlCompSequence(self, iCtrlComp):
    chemregDB, bcpvsDB = getDatabase(self)
    sSql = f"select pkey from {bcpvsDB}.CTRL_COMPOUND_SEQUENCE"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0]
    if int(iCtrlComp) > int(pkey):
        logger.info(f"Updating CTRL Comp counter to: {iCtrlComp}")
        sSql = f"update {bcpvsDB}.CTRL_COMPOUND_SEQUENCE set pkey={iCtrlComp}"
        cur.execute(sSql)

        
def isItNewStructure(self, molfile):
    """Return the compound_id of an existing match, '' if new, False on error."""
    chemregDB, bcpvsDB = getDatabase(self)
    std = standardizeMolecule(molfile)
    if not std.ok:
        return False
    mols = checkUniqueStructure(std.inchi_key, bcpvsDB)
    if not mols:
        return ''
    return mols[0][0]


def _getInchiKey(mol):
    """Helper to get InChI key for a mol object, or None on failure."""
    try:
        inchi = MolToInchi(mol, options='-w')
        if inchi:
            return InchiToInchiKey(inchi)
    except:
        pass
    return None


class StandardizedMolecule:
    """Single, deterministic result of running the registration
    standardization pipeline on an input molfile.

    Encapsulates everything downstream code needs so neither
    ChemRegAddMol nor BcpvsRegCompound has to call the RDKit tautomer
    canonicalizer twice on the same molecule (which used to be the
    source of duplicate compound rows because
    TautomerEnumerator.Canonicalize() is not deterministic across calls
    in RDKit 2022.3 for some scaffolds).
    """
    __slots__ = ('ok', 'std_molfile', 'std_smiles', 'inchi_key',
                 'tautomer_warning', 'error')

    def __init__(self, ok=False, std_molfile=None, std_smiles=None,
                 inchi_key=None, tautomer_warning='', error=''):
        self.ok = ok
        self.std_molfile = std_molfile
        self.std_smiles = std_smiles
        self.inchi_key = inchi_key
        self.tautomer_warning = tautomer_warning
        self.error = error


def standardizeMolecule(molfile, identifier=''):
    """Run the canonical standardization pipeline ONCE and return a
    `StandardizedMolecule` describing the result.

    Pipeline: parse -> Cleanup -> FragmentParent -> Uncharger ->
    TautomerEnumerator.Canonicalize -> molcart `mol2bin/bin2mol`
    roundtrip (so the molfile that goes into the structure tables
    matches what molcart will index).

    Two identity keys are derived:
      - `std_smiles`: canonical SMILES of the canonicalized tautomer
        (kept for the legacy `bcpvs.compound.smiles_std` column).
      - `inchi_key`: InChIKey computed AFTER cleanup/fragment/uncharge
        but BEFORE tautomer canonicalization. InChI normalizes
        tautomers itself, so the key is stable even when the RDKit
        tautomer canonicalizer flips. This is what registration uses
        to decide whether a compound is new.

    `tautomer_warning` is non-empty when the tautomer canonicalizer
    produced a different connectivity (same molecular formula) than
    the input. Registration is still allowed (this matches existing
    behavior) but the message is propagated to the client UI.
    """
    tag = f"[{identifier}] " if identifier else ''

    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return StandardizedMolecule(error=f'{tag}RDKit could not parse molfile')

    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False

    try:
        clean_mol = rdMolStandardize.Cleanup(mol, params)
    except Exception as e:
        logger.error(f"{tag}Cleanup failed: {e}")
        return StandardizedMolecule(error=f'{tag}Cleanup failed')
    try:
        parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    except Exception as e:
        logger.error(f"{tag}FragmentParent failed: {e}")
        return StandardizedMolecule(error=f'{tag}FragmentParent failed')

    uncharged_mol = rdMolStandardize.Uncharger().uncharge(parent_mol)

    # InChIKey BEFORE tautomer canonicalization. InChI is the stable
    # tautomer-aware identity used for de-duplication.
    pre_taut_key = _getInchiKey(uncharged_mol)

    try:
        taut_mol = rdMolStandardize.TautomerEnumerator(params).Canonicalize(uncharged_mol)
    except Exception as e:
        logger.error(f"{tag}TautomerEnumerator failed: {e}")
        return StandardizedMolecule(error=f'{tag}TautomerEnumerator failed')

    post_taut_key = _getInchiKey(taut_mol)

    tautomer_warning = ''
    if pre_taut_key and post_taut_key and pre_taut_key[:14] != post_taut_key[:14]:
        # Same input pre-tautomer InChIKey, different post. Same molecular
        # formula = tautomer flip; different formula = real connectivity
        # change which is suspicious and should also surface to the user.
        try:
            pre_mf = rdMolDescriptors.CalcMolFormula(uncharged_mol)
            post_mf = rdMolDescriptors.CalcMolFormula(taut_mol)
        except Exception:
            pre_mf = post_mf = ''
        if pre_mf == post_mf:
            tautomer_warning = (f'{tag}Tautomer canonicalization changed connectivity '
                                f'(same MF={pre_mf}); registered structure may differ '
                                f'from the supplier drawing.')
        else:
            tautomer_warning = (f'{tag}Tautomer canonicalization changed molecular '
                                f'formula ({pre_mf} -> {post_mf}); registered '
                                f'structure differs from the supplier drawing.')
        logger.warning(tautomer_warning)

    std_smiles = Chem.MolToSmiles(taut_mol)
    smiles_for_sql = std_smiles.replace('\\', '\\\\')

    # molcart roundtrip so the molfile we store matches what its
    # cartridge will index for substructure / similarity searches.
    cur.execute(f"select CONVERT(bin2mol(mol2bin(%s, 'smiles')) USING UTF8)",
                (smiles_for_sql,))
    saRes = cur.fetchall()
    if len(saRes) != 1 or saRes[0][0] is None:
        return StandardizedMolecule(error=f'{tag}molcart roundtrip failed')

    std_molfile = saRes[0][0]
    if isinstance(std_molfile, bytes):
        std_molfile = std_molfile.decode('utf-8', errors='replace')

    # Use the pre-tautomer InChIKey as the canonical identity, NOT
    # the post-tautomer one (the latter is what made registration
    # non-deterministic).
    inchi_key = pre_taut_key or post_taut_key

    return StandardizedMolecule(
        ok=True,
        std_molfile=std_molfile,
        std_smiles=std_smiles,
        inchi_key=inchi_key,
        tautomer_warning=tautomer_warning,
    )


def cleanStructureRDKit(molfile, identifier=''):
    """Backwards-compat shim that returns (std_molfile, std_smiles).

    New code should call `standardizeMolecule()` and use the returned
    `StandardizedMolecule.inchi_key` for compound lookups.
    """
    std = standardizeMolecule(molfile, identifier=identifier)
    if not std.ok:
        return False, False
    return std.std_molfile, std.std_smiles


def checkStructureIdentity(molfile, stdSMILES):
    """Compare original molfile with standardized SMILES to detect
    if tautomer canonicalization changed the molecular connectivity.
    Returns (is_match, original_smiles, details_string)"""
    try:
        orig_mol = Chem.MolFromMolBlock(molfile, removeHs=True, sanitize=True)
        if orig_mol is None:
            return True, '', 'Could not parse original molfile'

        params = rdMolStandardize.CleanupParameters()
        params.tautomerRemoveSp3Stereo = False
        params.tautomerRemoveBondStereo = False
        params.tautomerRemoveIsotopicHs = False
        try:
            clean_mol = rdMolStandardize.Cleanup(orig_mol, params)
            parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
            uncharger = rdMolStandardize.Uncharger()
            parent_mol = uncharger.uncharge(parent_mol)
        except:
            return True, '', 'Could not standardize original molfile for identity check'

        orig_smiles = Chem.MolToSmiles(parent_mol)

        std_mol = Chem.MolFromSmiles(stdSMILES)
        if std_mol is None:
            return True, orig_smiles, 'Could not parse stdSMILES'

        orig_inchi = MolToInchi(parent_mol, options='-w')
        std_inchi = MolToInchi(std_mol, options='-w')

        if orig_inchi is None or std_inchi is None:
            return True, orig_smiles, 'Could not generate InChI'

        orig_inchi_key = InchiToInchiKey(orig_inchi)
        std_inchi_key = InchiToInchiKey(std_inchi)

        if orig_inchi_key[:14] != std_inchi_key[:14]:
            # Check if this is just a tautomer (same molecular formula) or a real connectivity change
            orig_mf = rdMolDescriptors.CalcMolFormula(parent_mol)
            std_mf = rdMolDescriptors.CalcMolFormula(std_mol)

            if orig_mf == std_mf:
                # Same formula, different InChI = tautomer. Allow registration but log it.
                logger.info(f"checkStructureIdentity: tautomer difference allowed: "
                            f"original_smiles={orig_smiles}, std_smiles={stdSMILES}, "
                            f"MF={orig_mf}, orig_inchi_key={orig_inchi_key}, std_inchi_key={std_inchi_key}")
                return True, orig_smiles, ''

            details = (f"STRUCTURE MISMATCH: "
                      f"original_smiles={orig_smiles}, std_smiles={stdSMILES}, "
                      f"orig_mf={orig_mf}, std_mf={std_mf}, "
                      f"orig_inchi_key={orig_inchi_key}, std_inchi_key={std_inchi_key}")
            return False, orig_smiles, details

        return True, orig_smiles, ''
    except Exception as e:
        logger.error(f"checkStructureIdentity error: {str(e)}")
        return True, '', f'Error in structure identity check: {str(e)}'

def addStructure(database, molfile, newRegno, idColumnName):
    #####
    # Add the molecule to the structure tables
    cur.execute(
        f"INSERT INTO {database} (mol, {idColumnName}) VALUES (%s, %s)",
        (molfile, newRegno),
    )

    cur.execute(
        f"INSERT INTO {database}_MOL (mol, {idColumnName}) VALUES (mol2bin(%s, 'mol'), %s)",
        (molfile, newRegno),
    )

    cur.execute(
        f"""INSERT INTO {database}_ukey
            SELECT {idColumnName}, uniquekey(mol) AS molkey,
                   uniquekey(`mol`,'nostereo') AS molkeyns,
                   uniquekey(`mol`,'cistrans') AS molkeyct
            FROM {database} WHERE {idColumnName} = %s""",
        (newRegno,),
    )

    cur.execute(
        f"""INSERT INTO {database}_MOL_keysim
            SELECT {idColumnName}, fp(mol, 'sim') AS molkey
            FROM {database}_MOL WHERE {idColumnName} = %s""",
        (newRegno,),
    )

    cur.execute(
        f"""INSERT INTO {database}_MOL_key
            SELECT {idColumnName}, fp(mol, 'sss') AS molkey
            FROM {database}_MOL WHERE {idColumnName} = %s""",
        (newRegno,),
    )
    #
    #####

def registerNewCompound(bcpvsDB,
                        compound_id_numeric,
                        molfile,
                        stdSmiles,
                        inchiKey,
                        mf,
                        sep_mol_monoiso_mass,
                        ip_rights='',
                        compound_name=''):
    compound_id = f'CBK{compound_id_numeric}'
    sSql = f'''INSERT INTO {bcpvsDB}.compound (
        compound_id,
        compound_id_numeric,
        created_date,
        mf,
        ip_rights,
        sep_mol_monoiso_mass,
        smiles_std,
        smiles_std_string,
        inchi_key)
        VALUES (%s, %s, now(), %s, %s, %s, %s, %s, %s)'''
    cur.execute(sSql, (compound_id,
                       compound_id_numeric,
                       mf,
                       ip_rights,
                       sep_mol_monoiso_mass,
                       stdSmiles,
                       stdSmiles,
                       inchiKey))
    return compound_id


def registerNewBatch(bcpvsDB,
                     compound_id,
                     chemreg_regno,
                     notebook_ref,
                     suffix,
                     submitter,
                     project,
                     supplier,
                     biological_mw,
                     library_id,
                     compound_type,
                     product_type,
                     supplier_id,
                     supplier_batch,
                     restriction_comment='',
                     purity=-1):
    
    # --- Clean up default and empty values ---
    if not library_id:
        library_id = 'Lib-3002'
        
    # Convert empty suffix to Python's None, which MySQL translates to NULL
    if suffix == '':
        suffix = None

    # --- Parameterized SQL Query ---
    # Note: We still use an f-string for {bcpvsDB} because database/table names 
    # cannot be parameterized, only the actual data values can use %s placeholders.
    sSql = f'''INSERT INTO {bcpvsDB}.batch (
        compound_id,
        notebook_ref,
        suffix,
        submitter,
        submittal_date,
        project,
        purity,
        supplier,
        biological_mw,
        library_id,
        compound_type,
        product_type,
        supplier_id,
        supplier_batch,
        chemreg_regno,
        restriction_comment
    ) VALUES (
        %s, %s, %s, %s, now(), %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
    )'''
    
    # --- Tuple of values matching the %s placeholders in exact order ---
    values = (
        compound_id,
        notebook_ref,
        suffix,
        submitter,
        project,
        purity,
        supplier,
        biological_mw,
        library_id,
        compound_type,
        product_type,
        supplier_id,
        supplier_batch,
        chemreg_regno,
        restriction_comment
    )

    # --- Execute and handle errors ---
    try:
        cur.execute(sSql, values)
    except Exception as e:
        logger.error(f"Failed SQL: {sSql.strip()} | Values: {values}")
        logger.error(f"Error Details: {str(e)}")
        

def registerNewBatch_old(bcpvsDB,
                     compound_id,
                     chemreg_regno,
                     notebook_ref,
                     suffix,
                     submitter,
                     project,
                     supplier,
                     biological_mw,
                     library_id,
                     compound_type,
                     product_type,
                     supplier_id,
                     supplier_batch,
                     restriction_comment = '',
                     purity = -1):
    if not library_id:
        library_id = 'Lib-3002'
    sSql = f'''insert into {bcpvsDB}.batch (
    compound_id,
    notebook_ref,
    suffix,
    submitter,
    submittal_date,
    project,
    purity,
    supplier,
    biological_mw,
    library_id,
    compound_type,
    product_type,
    supplier_id,
    supplier_batch,
    chemreg_regno,
    restriction_comment)
    values (
    '{compound_id}',
    '{notebook_ref}',
    '{suffix}',
    '{submitter}',
    now(),
    '{project}',
    {purity},
    '{supplier}',
    {biological_mw},
    '{library_id}',
    '{compound_type}',
    '{product_type}',
    '{supplier_id}',
    '{supplier_batch}',
    {chemreg_regno},
    '{restriction_comment}')
    '''

    try:
        cur.execute(sSql)
    except Exception as e:
        logger.error(f"{sSql}")
        logger.error(f"{str(e)}")

    if suffix == '':
        sNULL = 'NULL'
        sSql = f"""update {bcpvsDB}.batch set
        suffix = {sNULL}
        where notebook_ref = '{notebook_ref}'"""
        cur.execute(sSql)


def getMoleculeProperties(self, molfile, chemregDB):
    sSql = f'''select
    bin2smiles(mol2bin('{molfile}'), 'mol') smiles,
    MolFormula(mol2bin('{molfile}', 'mol')),
    MolWeight(mol2bin(UNIQUEKEY('{molfile}', 'cistrans'))),
    MolNofMol(mol2bin('{molfile}', 'mol'))
    '''
    cur.execute(sSql)
    res = cur.fetchall()
    try:
        sSmiles = (res[0][0]).decode()
        C_MF = res[0][1].decode().replace(" ", "")
        mainFragMolWeight = res[0][2]
        iNrOfFragments = res[0][3]
    except:
        return (False, False, False, False, False, f'No smiles for molecule')

    try:
        molmassFormula = Formula(C_MF.replace('-', ''))
        C_CHNS = getAtomicComposition(molmassFormula.composition())
    except Exception as e:
        return (False, False, False, False, False, f'{str(e)}')

    if sSmiles == '':
        return (False, False, False, False, False, 'Empty molfile')
    if iNrOfFragments > 1:
        saSmileFragments = sSmiles.split('.')
    else:
        saSmileFragments = [sSmiles]

    if len(saSmileFragments) > 1 and iNrOfFragments == len(saSmileFragments):
        iFragPosition = 0
        for fragment in saSmileFragments:
            sSql = f'''select MolWeight(mol2bin(UNIQUEKEY('{fragment}', 'cistrans')))'''
            cur.execute(sSql)
            res = cur.fetchall()
            if res[0][0] == mainFragMolWeight:
                if iFragPosition != 0:
                    saSmileFragments[0], saSmileFragments[iFragPosition] = saSmileFragments[iFragPosition], saSmileFragments[0]
                    break
            iFragPosition += 1

    saSmileFragments = sorted(saSmileFragments, key=len, reverse=True)
    mainMolSmile = sSmiles
    saltSmile = ''

    if len(saSmileFragments) > 1:
        mainMolSmile = saSmileFragments[0]
        saltSmile = '.'.join(saSmileFragments[1:])

    cur.execute(f"""select
    MolWeight(mol2bin('{sSmiles}', 'smiles')),
    MolWeight(mol2bin('{mainMolSmile}', 'smiles'))
    """)
    resMolcart = cur.fetchall()
    C_MW = Formula(C_MF).mass
    C_MONOISO = resMolcart[0][1]
    if C_MONOISO == None:
        C_MONOISO = C_MW
    return (C_MF, C_MW, C_MONOISO, C_CHNS, saltSmile, '')


class PingDB(tornado.web.RequestHandler):
    def get(self):
        sSql = "select * from hive.user_details where pkey = 0"
        ret = cur.ping(sSql)
        if ret == 'error':
            self.set_status(400)

    def head(self):
        sSql = "select * from hive.user_details where pkey = 0"
        ret = cur.ping(sSql)
        if ret == 'error':
            self.set_status(400)


@jwtauth
class AddCtrlMol(tornado.web.RequestHandler):
    def post(self):
        compound_id = 'CBK999999'
        chemregDB, bcpvsDB = getDatabase(self)
        jpage = self.get_body_argument('jpage')
        chemist = self.get_body_argument('chemist')
        compound_type = self.get_body_argument('compound_type')
        project = self.get_body_argument('project')
        source = self.get_body_argument('source')
        solvent = self.get_body_argument('solvent')
        product = self.get_body_argument('product')
        library_id = self.get_body_argument('library_id')
        iCtrlComp =  self.get_body_argument('iCtrlComp')
        if updateCtrlCompSequence(self, iCtrlComp) == False:
            logger.error(f"Failed to update CTRL Comp counter")
            return
        sLib = re.search(r"(Lib-\d\d\d\d)", library_id)
        library_id = sLib.group()
        if library_id.startswith("Lib"):
            pass
        else:
            # If the library is blank set it to Compound collection lib
            library_id = 'Lib-3002'
        external_id = self.get_body_argument('external_id')
        supplier_batch = self.get_body_argument('supplier_batch')
        purity = self.get_body_argument('purity')
        ip_rights = self.get_body_argument('ip_rights')
        mw = self.get_body_argument('mw')
        restriction_comment = self.get_body_argument('restriction_comment')
        newRegno = getNewRegno(chemregDB)

        sSql = f"""
        insert into {chemregDB}.chem_info (
        regno,
        jpage,
        compound_id,
        rdate,
        chemist,
        compound_type,
        project,
        source,
        solvent,
        product,
        library_id,
        external_id,
        supplier_batch,
        purity,
        ip_rights,
        C_MW)
        values (
        '{newRegno}',
        '{jpage}',
        '{compound_id}',
        now(),
        '{chemist}',
        '{compound_type}',
        '{project}',
        '{source}',
        '{solvent}',
        '{product}',
        '{library_id}',
        '{external_id}',
        '{supplier_batch}',
        {purity},
        '{ip_rights}',
        {mw}
        )
        """

        try:
            cur.execute(sSql)
            pass
        except Exception as e:
            logger.error(f"{sSql}")
            logger.error(f"{str(e)}")
        
        saSalts = ''
        registerNewBatch(bcpvsDB,
                         compound_id,
                         newRegno,
                         jpage,
                         saSalts,
                         chemist,
                         project,
                         source,         #  Supplier
                         mw,
                         library_id,
                         compound_type,
                         product,
                         external_id,
                         supplier_batch,
                         restriction_comment,
                         purity = -1)


@jwtauth
class AddNostructMol(tornado.web.RequestHandler):
    def post(self):
        compound_id = 'CBK999999'
        chemregDB, bcpvsDB = getDatabase(self)
        jpage = self.get_body_argument('jpage')
        chemist = self.get_body_argument('chemist')
        compound_type = self.get_body_argument('compound_type')
        project = self.get_body_argument('project')
        source = self.get_body_argument('source')
        solvent = self.get_body_argument('solvent')
        product = self.get_body_argument('product')
        library_id = self.get_body_argument('library_id')
        sLib = re.search(r"(Lib-\d\d\d\d)", library_id)
        library_id = sLib.group()
        if library_id.startswith("Lib"):
            pass
        else:
            # If the library is blank set it to Compound collection lib
            library_id = 'Lib-3002'
        external_id = self.get_body_argument('external_id')
        supplier_batch = self.get_body_argument('supplier_batch')
        purity = self.get_body_argument('purity')
        ip_rights = self.get_body_argument('ip_rights')
        mw = self.get_body_argument('mw')
        restriction_comment = self.get_body_argument('restriction_comment')
        newRegno = getNewRegno(chemregDB)

        
        sSql = f"""
        insert into {chemregDB}.chem_info (
        regno,
        jpage,
        compound_id,
        rdate,
        chemist,
        compound_type,
        project,
        source,
        solvent,
        product,
        library_id,
        external_id,
        supplier_batch,
        purity,
        ip_rights,
        C_MW)
        values (
        '{newRegno}',
        '{jpage}',
        '{compound_id}',
        now(),
        '{chemist}',
        '{compound_type}',
        '{project}',
        '{source}',
        '{solvent}',
        '{product}',
        '{library_id}',
        '{external_id}',
        '{supplier_batch}',
        {purity},
        '{ip_rights}',
        {mw}
        )
        """

        try:
            cur.execute(sSql)
            pass
        except Exception as e:
            logger.error(f"{sSql}")
            logger.error(f"{str(e)}")
        
        saSalts = ''
        registerNewBatch(bcpvsDB,
                         compound_id,
                         newRegno,
                         jpage,
                         saSalts,
                         chemist,
                         project,
                         source,         #  Supplier
                         mw,
                         library_id,
                         compound_type,
                         product,
                         external_id,
                         supplier_batch,
                         restriction_comment,
                         purity = -1)


@jwtauth
class GetMolkeyct(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f'''
        select CONVERT(molkeyct USING utf8) from {bcpvsDB}.JCMOL_MOLTABLE_ukey
        where compound_id = '{sCmpId}'
        '''
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            res = res[0][0]
        else:
            res = ''
        self.finish(res)


@jwtauth
class GetCompoundDuplicates(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f'''
        select compound_id from {bcpvsDB}.JCMOL_MOLTABLE_ukey
        where molkeyct = (select molkeyct from {bcpvsDB}.JCMOL_MOLTABLE_ukey
                         where compound_id = '{sCmpId}')
        '''
        cur.execute(sSql)
        res = res2json()
        #res = cur.fetchall()
        if len(res) > 0:
            pass
            #res = res2json()
        else:
            res = ''
        self.finish(res)


@jwtauth
class CreateSalt(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sSmiles = self.get_argument("smiles")
        iSaltNumber = getNewSaltNumber(chemregDB)
        mol = Chem.MolFromSmiles(sSmiles)
        mw = Descriptors.MolWt(mol)
        mf = rdMolDescriptors.CalcMolFormula(mol)
        suffix = f'X{iSaltNumber}'
        sio = StringIO()
        with Chem.SDWriter(sio) as w:
            w.write(mol)
        molfile = sio.getvalue()
        sSql = f"""insert into {chemregDB}.salts (suffix, smiles, mw, mf, molfile)
        values ('{suffix}', '{sSmiles}', {mw}, '{mf}', '{molfile}')
        """
        cur.execute(sSql)


@jwtauth
class GetCanonicSmiles(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sSmiles = self.get_argument("smiles")
        sLetters, saRemainderSmile = getSaltLetters([sSmiles], chemregDB)
        mol = Chem.MolFromSmiles(sSmiles)
        flattenSmile = Chem.rdmolfiles.MolToSmiles(mol)
        canonSmile = Chem.CanonSmiles(flattenSmile)
        sRes = json.dumps([f'{sLetters}', canonSmile])
        self.finish(sRes)


@jwtauth
class GetLastCtrlComp(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sEln = self.get_argument("eln")
        eln_id = sEln[0:6]
        sSql = f'''select pkey from {bcpvsDB}.CTRL_COMPOUND_SEQUENCE'''
        cur.execute(sSql)
        cRes = cur.fetchall()
        iBatch = int(cRes[0][0])
        self.finish(f'{iBatch}')


@jwtauth
class GetLastBatchFromEln(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sEln = self.get_argument("eln")
        eln_id = sEln[0:6]
        sSql = f'''select notebook_ref from {bcpvsDB}.batch
        where notebook_ref like '{eln_id}%' order by notebook_ref desc'''
        cur.execute(sSql)
        cRes = cur.fetchall()
        if len(cRes) > 0:
            iBatch = int((cRes[0][0])[-3:])
            self.finish(f'{iBatch}')
        else:
            self.finish(b'0')


@jwtauth
class GetNextSdfSequence(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sSql = f"select pkey from {chemregDB}.sdfile_sequence"
        cur.execute(sSql)
        pkey = cur.fetchall()[0][0] +1
        sSql = f"update {chemregDB}.sdfile_sequence set pkey={pkey}"
        cur.execute(sSql)
        self.finish(f'{pkey}')


@jwtauth
class UpdateStructureAdmin(tornado.web.RequestHandler):
    def post(self):
        chemregDB, bcpvsDB = getDatabase(self)
        molfile = self.get_body_argument('molfile')
        compound_id = self.get_body_argument('compound_id')
        smiles = self.get_body_argument('smiles')
        database = bcpvsDB + '.JCMOL_MOLTABLE'
        (C_MF,
         C_MW,
         C_MONOISO,
         C_CHNS,
         saSalts,
         errorMessage) = getMoleculeProperties(self, molfile, chemregDB)
        sSql = f'''
        update {bcpvsDB}.compound
        set
        mf = '{C_MF}',
        sep_mol_monoiso_mass = '{C_MONOISO}',
        smiles_std = '{smiles}',
        smiles_std_string = '{smiles}'
        where compound_id = '{compound_id}'
        '''
        cur.execute(sSql)

        sSql = f'''
update {database}
set mol = '{molfile}'
where compound_id = '{compound_id}'
        '''
        cur.execute(sSql)
        ######################################################
        ###  Update the molcart tables for compound

        sSql = f"""
update {database}_MOL
set mol = mol2bin('{molfile}', 'mol')
where compound_id = '{compound_id}'
        """
        cur.execute(sSql)

        sSql = f"""
update {database}_ukey
set
molkey = (select uniquekey(mol)
from {database} where compound_id = '{compound_id}'),
molkeyns = (select uniquekey(`mol`,'nostereo')
from {database} where compound_id = '{compound_id}'),
molkeyct = (select uniquekey(`mol`,'cistrans')
from {database} where compound_id = '{compound_id}')
where compound_id = '{compound_id}'
        """
        cur.execute(sSql)

        sSql = f"""
update {database}_MOL_keysim set molkey = (select fp(mol, 'sim')
from {database}_MOL where compound_id = '{compound_id}')
where compound_id = '{compound_id}'
        """
        cur.execute(sSql)

        sSql = f"""
update {database}_MOL_key set molkey = (select fp(mol, 'sss') as molkey
from {database}_MOL where compound_id = '{compound_id}')
where compound_id = '{compound_id}'
        """
        cur.execute(sSql)



@jwtauth
class GetRegnosFromSdfSequence(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sdfile_sequence = self.get_argument("sdfile_sequence")
        sSql = f'''select regno from {chemregDB}.chem_info
        where sdfile_sequence = {sdfile_sequence}
        '''
        cur.execute(sSql)
        res = res2json()
        self.finish(res)


@jwtauth
class BcpvsRegCompound(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        cur.execute(f'''SELECT regno, compound_id, c_mf, ip_rights, c_monoiso,
                              molfile, jpage, suffix, chemist, project, source,
                              c_mw, library_id, compound_type, product,
                              external_id, supplier_batch, purity, inchi_key
                       FROM {chemregDB}.chem_info WHERE regno = %s''',
                    (regno,))
        row = cur.fetchall()
        if not row:
            self.set_status(404)
            self.finish(f'No chem_info row for regno {regno}'.encode())
            return
        (regno,
         compound_id,
         c_mf,
         ip_rights,
         c_monoiso,
         molfile,
         jpage,
         suffix,
         chemist,
         project,
         source,
         c_mw,
         library_id,
         compound_type,
         product,
         external_id,
         supplier_batch,
         purity,
         inchi_key) = row[0]

        # Idempotency: if a batch already exists for this regno + jpage,
        # do not create a duplicate batch row. (The retry path in the
        # client iterates regnos one-by-one and a transient network
        # error used to register the same batch twice.)
        cur.execute(f'''SELECT compound_id FROM {bcpvsDB}.batch
                       WHERE chemreg_regno = %s AND notebook_ref = %s''',
                    (regno, jpage))
        existing = cur.fetchall()
        if existing:
            existing_compound_id = existing[0][0]
            logger.info(f"BcpvsRegCompound: batch for regno={regno} jpage={jpage} "
                        f"already exists with compound_id={existing_compound_id}, "
                        f"skipping duplicate registration.")
            self.finish(b'compound_id')
            return

        if not inchi_key:
            # Older row (pre-migration or backfill not yet run). Compute
            # standardization once here and persist the inchi_key for
            # future calls.
            std = standardizeMolecule(molfile, identifier=f'regno={regno} jpage={jpage}')
            if not std.ok:
                self.set_status(500)
                self.finish(f'RDKit standardization failed for regno {regno}: {std.error}'.encode())
                return
            inchi_key = std.inchi_key
            std_smiles = std.std_smiles
            std_molfile = std.std_molfile
            cur.execute(f'UPDATE {chemregDB}.chem_info SET inchi_key = %s WHERE regno = %s',
                        (inchi_key, regno))
        else:
            # ChemRegAddMol already standardized; we trust the persisted
            # inchi_key. We still need the std_smiles + std_molfile for
            # the bcpvs.compound row if we end up creating it.
            std = standardizeMolecule(molfile, identifier=f'regno={regno} jpage={jpage}')
            if not std.ok:
                self.set_status(500)
                self.finish(f'RDKit standardization failed for regno {regno}: {std.error}'.encode())
                return
            std_smiles = std.std_smiles
            std_molfile = std.std_molfile

        (C_MF, C_MW, C_MONOISO, C_CHNS, saSalts, errorMessage) = \
            getMoleculeProperties(self, molfile, chemregDB)

        if compound_id in ('', None):
            mols = checkUniqueStructure(inchi_key, bcpvsDB)
            if mols:
                compound_id = mols[0][0]
            else:
                # Register a new compound
                compound_id_numeric = getNewCompoundId(bcpvsDB)
                compound_id = registerNewCompound(bcpvsDB,
                                                  compound_id_numeric,
                                                  std_molfile,
                                                  std_smiles,
                                                  inchi_key,
                                                  C_MF,
                                                  C_MONOISO,
                                                  ip_rights,
                                                  compound_name='')
                addStructure(f"{bcpvsDB}.JCMOL_MOLTABLE",
                             std_molfile,
                             compound_id,
                             'compound_id')

        cur.execute(f'UPDATE {chemregDB}.chem_info SET compound_id = %s WHERE regno = %s',
                    (compound_id, regno))

        registerNewBatch(bcpvsDB,
                         compound_id,
                         regno,
                         jpage,
                         saSalts,
                         chemist,
                         project,
                         source,
                         c_mw,
                         library_id,
                         compound_type,
                         product,
                         external_id,
                         supplier_batch,
                         restriction_comment='',
                         purity=-1)
        self.finish(b'compound_id')


@jwtauth
class ChemRegAddMol(tornado.web.RequestHandler):
    def post(self):
        chemregDB, bcpvsDB = getDatabase(self)
        molfile = self.get_body_argument('molfile')
        jpage = self.get_body_argument('jpage')
        chemist = self.get_body_argument('chemist')
        compound_type = self.get_body_argument('compound_type')
        project = self.get_body_argument('project')
        source = self.get_body_argument('source')
        solvent = self.get_body_argument('solvent')
        product = self.get_body_argument('product')
        library_id = self.get_body_argument('library_id')
        sLib = re.search(r"(Lib-\d\d\d\d)", library_id)
        # library_id is required by the schema but the user does not
        # always supply one; fall back to the catch-all "Compound
        # collection" library Lib-3002.
        library_id = sLib.group() if sLib else 'Lib-3002'

        jpage = re.sub(r"\s+", "", jpage)

        external_id = self.get_body_argument('external_id')
        supplier_batch = self.get_body_argument('supplier_batch')
        purity = self.get_body_argument('purity')
        ip_rights = self.get_body_argument('ip_rights')
        sdfile_sequence = self.get_body_argument('sdfile_sequence')

        if "0  0  0     0  0            999 V2" in molfile:
            self.set_status(500)
            self.finish(f'Nostruct for jpage: {jpage}')
            return

        # ---------- Standardize ONCE ----------
        std = standardizeMolecule(molfile, identifier=f'jpage={jpage}')
        if not std.ok:
            self.set_status(500)
            self.finish(f'RDKit failed to convert Molfile for: {jpage} ({std.error})')
            return

        # ---------- Validate BEFORE any INSERT ----------
        is_match, orig_smiles, mismatch_details = checkStructureIdentity(molfile, std.std_smiles)
        if not is_match:
            logger.error(f"ChemRegAddMol: Structure identity FAILED for jpage {jpage}: "
                         f"{mismatch_details}")
            self.set_status(500)
            self.finish(f'Structure mismatch for {jpage}: original={orig_smiles}, '
                        f'standardized={std.std_smiles}')
            return

        (C_MF, C_MW, C_MONOISO, C_CHNS, saSalts, errorMessage) = \
            getMoleculeProperties(self, molfile, chemregDB)
        if saSalts is None:
            saSalts = ''
        if C_MF is False:
            self.set_status(500)
            self.finish(f'Molfile failed {external_id} {errorMessage}')
            logger.error(f'Molfile failed {external_id} {errorMessage}')
            return
        if saSalts is False:
            self.set_status(500)
            self.finish(f'Unknown salt in molfile for {external_id} {errorMessage}')
            logger.error(f'Unknown salt in molfile for {external_id} {errorMessage}')
            return

        # ---------- Compound lookup by InChIKey ----------
        mols = checkUniqueStructure(std.inchi_key, bcpvsDB)
        if mols:
            sStatus = 'oldMolecule'
            compound_id = mols[0][0]
        else:
            sStatus = 'newMolecule'
            compound_id = ''

        if purity == '':
            purity = -1

        newRegno = getNewRegno(chemregDB)

        # NOTE: chem_info.MOLFILE intentionally stores the SUPPLIER's
        # original molfile (not std.std_molfile). The standardized form
        # used for matching is captured by chem_info.inchi_key.
        cur.execute(f"""INSERT INTO {chemregDB}.chem_info (
            regno, jpage, compound_id, rdate, chemist, compound_type, project,
            source, solvent, product, library_id, external_id, supplier_batch,
            purity, ip_rights, C_CHNS, C_MF, C_MW, C_MONOISO, SUFFIX,
            sdfile_sequence, molfile, inchi_key)
            VALUES (%s, %s, %s, now(), %s, %s, %s, %s, %s, %s, %s, %s, %s,
                    %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)""",
                    (newRegno, jpage, compound_id, chemist, compound_type, project,
                     source, solvent, product, library_id, external_id, supplier_batch,
                     purity, ip_rights, C_CHNS, C_MF, C_MW, C_MONOISO, saSalts,
                     sdfile_sequence, molfile, std.inchi_key))

        if saSalts == '':
            cur.execute(f"UPDATE {chemregDB}.chem_info SET suffix = NULL WHERE regno = %s",
                        (newRegno,))

        # CHEM structure tables get the standardized molfile so molcart
        # indexing is consistent across re-registrations.
        addStructure(f"{chemregDB}.CHEM", std.std_molfile, newRegno, 'regno')

        # Surface tautomer warnings to the client UI without failing the
        # registration (registration succeeded, but the user should
        # review the structure).
        if std.tautomer_warning:
            self.set_header('X-Tautomer-Warning', std.tautomer_warning[:512])
            self.finish(f'{sStatus};warning={std.tautomer_warning}')
            return

        self.finish(sStatus)


@jwtauth
class Search(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        value = self.get_argument("value")
        values = (value, )

        if column == 'JPAGE':
            sSql = f"select regno from {chemregDB}.chem_info where {column} like '{value}%'"
        else:
            sSql = f"select regno from {chemregDB}.chem_info where {column} = '{value}'"
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetRegnoData(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        regno = self.get_argument("regno")
        values = (regno, )
        sSql = f"select {column} from {chemregDB}.chem_info where regno = {regno}"
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class CreateMolImageFromMolfile(tornado.web.RequestHandler):
    def post(self):
        sStatus = ''
        molfile = '\n' + self.get_body_argument('molfile')
        sMolId = self.get_body_argument('mol_id')
        try:
            os.remove(f'mols/{sMolId}.png')
        except:
            pass
        if len(molfile) < 85:
            # The molecule is a nostruct, rdkit doesn't like that
            resDict = {
                "molfile": '',
                "status": '',
                "smiles": ''
            }
            self.finish(json.dumps(resDict))
            return
        mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
        params = rdMolStandardize.CleanupParameters()
        params.tautomerRemoveSp3Stereo = False
        params.tautomerRemoveBondStereo = False
        params.tautomerRemoveIsotopicHs = False
        try:
            clean_mol = rdMolStandardize.Cleanup(mol, params)
        except:
            logger.error(f'Failed to standardize molfile {sMolId}')
            sStatus = 'Molecule altered, checkit'
            sSql = f'''select bin2mol(moldepict(mol2bin(UNIQUEKEY('{molfile}',
                                                                  'cistrans'),
                                                                  'smiles')))'''
            cur.execute(sSql)
            molfile = cur.fetchall()[0][0].decode("utf-8")
            mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
            if mol == None:
                logger.error(f'Failed to standardize molfile: {molfile}')
                return
            params = rdMolStandardize.CleanupParameters()
            params.tautomerRemoveSp3Stereo = False
            params.tautomerRemoveBondStereo = False
            params.tautomerRemoveIsotopicHs = False
            clean_mol = rdMolStandardize.Cleanup(mol, params)

        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
        uncharger = rdMolStandardize.Uncharger()
        uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
        te = rdMolStandardize.TautomerEnumerator(params) # idem
        taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)
        smiles = Chem.MolToSmiles(taut_uncharged_parent_clean_mol)
        m = Chem.MolToMolBlock(taut_uncharged_parent_clean_mol)

        # Don't think we need to generate new coords with molcart
        #sSql = f'''select bin2mol(moldepict(mol2bin(UNIQUEKEY('{m}', 'cistrans'), 'smiles')))'''
        #cur.execute(sSql)
        #strip = cur.fetchall()[0][0].decode("utf-8")
        #m = Chem.MolFromMolBlock(strip)

        m = 'Id' + m
        strip = m
        new_mol = Chem.MolFromMolBlock(strip, removeHs=False, sanitize=True)
        try:
            Draw.MolToFile(new_mol, f'mols/{sMolId}.png', kekulize=True, size=(280, 280))
        except Exception as e:
            logger.error(str(e))
            logger.error(f"{sMolId} is nostruct, png failed")
        resDict = {
            "molfile": strip,
            "status": sStatus,
            "smiles": smiles
        }

        self.finish(json.dumps(resDict))

@jwtauth
class CreateMolImage(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""select molfile from {chemregDB}.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        molfile = cur.fetchall()

        if len(molfile) == 0:
            sSql = f"""select molfile from chemspec.chem_info
            where regno = '{regno}'"""
            cur.execute(sSql)
            molfile = cur.fetchall()

        if len(molfile) > 0 and molfile[0][0] != None:
            try:
                createPngFromMolfile(regno, molfile[0][0])
            except:
                pass
        self.finish()


@jwtauth
class CreateBcpvsMolImage(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f"""select mol, compound_id from {bcpvsDB}.JCMOL_MOLTABLE
                   where compound_id like '{sCmpId}%'"""
        cur.execute(sSql)
        molfile = cur.fetchall()
        if len(molfile) > 0 and molfile[0][0] != None:
            try:
                createPngFromMolfile(sCmpId, molfile[0][0])
            except:
                pass
        else:
            try:
                os.remove(f'mols/{sCmpId}.png')
            except:
                pass
        self.finish()


@jwtauth
class LoadMolfile(tornado.web.RequestHandler):
    def post(self):
        chemregDB, bcpvsDB = getDatabase(self)
        fBody = self.request.files['file'][0]
        regnoBody = self.request.files['regno'][0]
        molfile = tornado.escape.xhtml_unescape(fBody.body)
        regno = tornado.escape.xhtml_unescape(regnoBody.body)

        (C_MF,
         C_MW,
         C_MONOISO,
         C_CHNS,
         saSalts,
         errorMessage) = getMoleculeProperties(self, molfile, chemregDB)
        if C_MF == False:
            self.set_status(500)
            self.finish(f'Molfile failed {regno} {errorMessage}')
            return
        try:
            addStructure(f'{chemregDB}.CHEM', molfile, regno, 'regno')
        except Exception as e:
            logger.error(str(e))
            logger.error(f'failed on regno {regno}')
            return
        compound_id = isItNewStructure(self, molfile)
        sSql = f"""update {chemregDB}.chem_info set
                   `molfile` = '{molfile}',
                    compound_id = '{compound_id}',
                    C_MF = '{C_MF}',
                    C_MW = {C_MW},
                    C_MONOISO = {C_MONOISO},
                    C_CHNS = '{C_CHNS}',
                    suffix = '{saSalts}'
                  where regno = {regno}"""
        cur.execute(sSql)

        if saSalts == '':
            sNULL = 'NULL'
            sSql = f"""update {chemregDB}.chem_info set
            suffix = {sNULL}
            where regno = {regno}"""
            cur.execute(sSql)

        createPngFromMolfile(regno, molfile)


@jwtauth
class UpdateRegnoBatch(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        batch = self.get_argument("batch")
        sSql = f"""select compound_id from {bcpvsDB}.batch
                   where notebook_ref = '{batch}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.set_status(500)
            self.finish()
        else:
            sSql = f"""select regno from {chemregDB}.chem_info
            where jpage = '{batch}' and regno != '{regno}'"""
            cur.execute(sSql)
            res = cur.fetchall()
            if len(res) > 0:
                self.set_status(500)
                self.finish()
            else:
                sSql = f"""update {chemregDB}.chem_info
                set jpage = '{batch}' where regno = '{regno}'"""
                cur.execute(sSql)
                self.finish()


@jwtauth
class CreateSupplier(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        supplierName = self.get_argument("supplier")
        sSql = f"""insert into {bcpvsDB}.compound_suppliers
        (name) values ('{supplierName}')"""
        cur.execute(sSql)


@jwtauth
class CreateRegno(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""insert into {chemregDB}.chem_info
        (regno, rdate) values (%s, now())"""
        val = (regno, )
        cur.execute(sSql, val)


@jwtauth
class DeleteRegno(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""select regno, c_mf, chemist from {chemregDB}.chem_info
                  where c_mf is null and chemist is null and
                  regno={regno}"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            logger.info('Deleting ' + str(res))
            sSql = f"""delete from {chemregDB}.chem_info
                      where regno = {regno}"""
            cur.execute(sSql)


@jwtauth
class UpdateColumn(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        value = self.get_argument("value")
        regno = self.get_argument("regno")

        if 'PURITY' in column:
            try:
                value = int(value)
            except:
                logger.info(f"Purity is non integer, setting purity to -1")
                value = -1
        values = (value, regno, )
        sSql = f"""update {chemregDB}.chem_info set {column} = '{value}'
        where regno = {regno}"""

        if column == 'library_id':
            sLib = re.search(r"(Lib-\d\d\d\d)", value)
            try:
                value = sLib.group()
                sSql = f"""update {chemregDB}.chem_info set {column} = '{value}'
                where regno = {regno}"""
            except:
                sSql = f"""update {chemregDB}.chem_info set {column} = NULL
                where regno = {regno}"""

        cur.execute(sSql)


@jwtauth
class RegCompound(tornado.web.RequestHandler):
    def get(self):
        # Contains user found in previous auth
        if self.request.headers.get('auth'):
            self.write('ok')


@jwtauth
class GetCompound(tornado.web.RequestHandler):
    def get(self):
        # Contains user found in previous auth
        if self.request.headers.get('auth'):
            self.write('ok')


@jwtauth
class GetMolfile(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""select molfile from {chemregDB}.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(res[0][0])
            return

        sSql = f"""select molfile from chemspec.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(res[0][0])


@jwtauth
class GetForwardCompound(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f"""select compound_id from {bcpvsDB}.compound
                   where compound_id > '{sCmpId}' order by compound_id asc limit 1"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetBackwardCompound(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f"""select compound_id from {bcpvsDB}.compound
                   where compound_id < '{sCmpId}' order by compound_id desc limit 1"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetBackwardRegno(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")

        sSql = f"""select greatest(chemreg, chemspec) from (
        WITH
        cte1 AS (SELECT IFNULL((SELECT regno FROM {chemregDB}.chem_info
        where regno < '{regno}' order by regno desc limit 1), 0) as chemreg),
        cte2 AS (select IFNULL((SELECT regno chemspec FROM chemspec.chem_info
        where regno < '{regno}' order by regno desc limit 1), 0) as chemspec)
        SELECT chemreg, chemspec FROM cte1, cte2) ss
        """
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetForwardRegno(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")

        sSql = f"""select least(chemreg, chemspec) from (
        WITH
        cte1 AS (SELECT regno chemreg FROM {chemregDB}.chem_info
        where regno > '{regno}' order by regno asc limit 1),
        cte2 AS (select IFNULL((SELECT regno chemspec FROM chemspec.chem_info
        where regno > '{regno}' order by regno asc limit 1), 999999999) as chemspec)
        SELECT chemreg, chemspec FROM cte1, cte2) ss
        """
        cur.execute(sSql)
        res = cur.fetchall()

        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetRegnoFromCompound(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        if sCmpId.startswith('CBK6'):
            sSql = f"""select regno from {chemregDB}.chem_info
                       where compound_id = '{sCmpId}'"""
        else:
            sSql = f"""
            select cs.regno from chemspec.chem_info cs
            join {bcpvsDB}.batch b on cs.regno = b.chemspec_regno
            where b.compound_id = '{sCmpId}'
                    """

        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetCompoundFromRegno(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sRegno = self.get_argument("regno")
        sSql = f"""select compound_id from {chemregDB}.chem_info
                   where regno = '{sRegno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))
            return

        sSql = f"""
        select compound_id from {bcpvsDB}.batch b
        join chemspec.chem_info ci on ci.regno = b.chemspec_regno
        where regno = '{sRegno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetMolfileBcpvs(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        compound_id = self.get_argument("compound_id")
        sSql = f"""select mol from {bcpvsDB}.JCMOL_MOLTABLE
                   where compound_id = '{compound_id}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(res[0][0])


@jwtauth
class GetTextColumn(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        regno = self.get_argument("regno")
        if regno == None:
            return

        if column == 'library_id':
            sSql = f"""
select IFNULL((select concat(library_id, " ", l.description)
from {chemregDB}.chem_info c, {bcpvsDB}.compound_library l
where c.library_id = l.library_name
and c.regno = '{regno}'), null)
            """
        else:
            sSql = f"""select cast({column} as char)
            from {chemregDB}.chem_info where regno = {regno}"""
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetColComboData(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        sSql = ''
        if column == 'project':
            sSql = """select project_name from hive.project_details
                      order by project_name"""
            sSql = """select project_name from hive.project_details
                      order by created_date desc"""
        elif column == 'chemist':
            sSql = """select fullname from hive.user_details
            where ORGANIZATION = 'chemistry'"""
        elif column == 'compound_type':
            sSql = f"select type from {bcpvsDB}.compound_type order by type"
        elif column == 'product':
            sSql = f"SELECT type FROM {bcpvsDB}.product_type order by type"
        elif column == 'supplier':
            sSql = f"SELECT name FROM {bcpvsDB}.compound_suppliers order by name"
        elif column == 'solvent':
            sSql = "SELECT solvent FROM chemspec.solvent_tbl order by solvent"
        elif column == 'library_id':
            sSql = f"""SELECT concat(library_name, " ", description) as library_name
                       FROM {bcpvsDB}.compound_library order by library_name"""
        elif column == 'library_description':
            sSql = f"""SELECT description FROM {bcpvsDB}.compound_library
                      order by description"""
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetLibraryName(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        try:
            library_id = self.get_argument("library_id")
        except:
            library_id = ''
        sLib = re.search(r"(Lib-\d\d\d\d)", library_id)
        if sLib == None:
            self.write(b'[]')
            return
        library_id = sLib.group()
        sSql = f"""SELECT description FROM {bcpvsDB}.compound_library
                  where library_name = '{library_id}'"""
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetLibraries(tornado.web.RequestHandler):
    def get(self):
        sSql = "select fullname from hive.user_details where ORGANIZATION = 'chemistry'"
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetNextRegno(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        newRegno = getNewRegno(chemregDB)
        self.write(json.dumps(newRegno))


@jwtauth
class CreateLibrary(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sLibraryDescription = self.get_argument("library_name")
        sSupplier = self.get_argument("supplier")
        restricted = self.get_argument("restricted")
        if restricted == 'False':
            restricted = 'No'
        else:
            restricted = 'Yes'
            
        sSql = f"select pkey from {bcpvsDB}.library_id_sequence"
        cur.execute(sSql)
        pkey = cur.fetchall()[0][0] +1
        sSql = f"update {bcpvsDB}.library_id_sequence set pkey={pkey}"
        cur.execute(sSql)

        library_id = f'Lib-{pkey}'

        sSql = f"""insert into {bcpvsDB}.compound_library
        (library_name,
        supplier,
        description,
        restricted) values
        ('{library_id}',
        '{sSupplier}',
        '{sLibraryDescription}',
        '{restricted}')
        """
        cur.execute(sSql)


class GetDatabase(tornado.web.RequestHandler):
    def get(self):
        sRes = json.dumps([['Live'], ['Test']])
        self.finish(sRes)

@jwtauth
class UploadBinary(tornado.web.RequestHandler):
    def post(self, *args, **kwargs):
        os_name = self.get_argument("os_name")

        try:
            # self.request.files['file'][0]:
            # {'body': 'Label Automator ___', 'content_type': u'text/plain', 'filename': u'k.txt'}
            file1 = self.request.files['file'][0]
        except:
            logging.error("Error cant find file1 in the argument list")
            return

        bin_file = ""
        if os_name == 'Windows':
            bin_file = f'dist/{os_name}/ch.exe'
        elif os_name == 'Linux':
            bin_file = f'dist/{os_name}/ch'
        elif os_name == 'Darwin':
            bin_file = f'dist/{os_name}/ch'
        else:
            # unsupported OS
            self.set_status(500)
            self.write({'message': 'OS not supported'})
            return

        output_file = open(bin_file, 'wb')
        output_file.write(file1['body'])
        output_file.close()

class GetChemRegBinary(tornado.web.RequestHandler):
    def post(self):
        pass

    def get(self, os_name):
        bin_file = ""
        if os_name == 'Windows':
            bin_file = f'dist/{os_name}/ch.exe'
        elif os_name == 'Linux':
            bin_file = f'dist/{os_name}/ch'
        elif os_name == 'Darwin':
            bin_file = f'dist/{os_name}/ch'
        else:
            # unsupported OS
            self.set_status(500)
            self.write({'message': 'OS not supported'})
            return
        try:
            with open(bin_file, 'rb') as f:
                logging.info("sending bin file")
                self.set_status(200)
                self.write(f.read())
        except Exception as e:
            logging.error(f"Did not send bin file, error: {str(e)}")

class getVersionData(tornado.web.RequestHandler):
    def post(self):
        pass

    def get(self):
        try:
            with open('./ver.dat', 'r') as f:
                self.write(json.load(f))
                return
        except Exception as e:
            logging.error(str(e))
            self.set_status(500)
            self.write({'message': 'ver.dat not available'})

@jwtauth
class UploadVersionNo(tornado.web.RequestHandler):
    def post(self, *args, **kwargs):
        ver_no = self.get_argument("ver_no")
        ver_file = "ver.dat"

        with open(ver_file, "r") as f:
            data = json.load(f)

        data["version"] = ver_no

        with open(ver_file, "w") as f:
            json.dump(data, f)

@jwtauth
class UploadLauncher(tornado.web.RequestHandler):
    def post(self, *args, **kwargs):
        os_name = self.get_argument("os_name")

        try:
            # self.request.files['file'][0]:
            # {'body': 'Label Automator ___', 'content_type': u'text/plain', 'filename': u'k.txt'}
            file1 = self.request.files['file'][0]
        except:
            logging.error("Error cant find file1 in the argument list")
            return

        bin_file = ""
        if os_name == 'Windows':
            bin_file = f'dist/launchers/{os_name}/chemreg.exe'
        elif os_name == 'Linux':
            bin_file = f'dist/launchers/{os_name}/chemreg'
        elif os_name == 'Darwin':
            bin_file = f'dist/launchers/{os_name}/chemreg'
        else:
            # unsupported OS
            self.set_status(500)
            self.write({'message': 'OS not supported'})
            return

        output_file = open(bin_file, 'wb')
        output_file.write(file1['body'])
        output_file.close()


class ChemblExport(tornado.web.RequestHandler):
    def get(self, sRIDX, sBatches):
        pass
    def post(self, *args, **kwargs):
        try:
            sRIDX = self.get_argument("RIDX").strip()
            sAIDX = self.get_argument("AIDX").strip()
            sProject = self.get_argument("project").strip()
            sELN = self.get_argument("ELN").strip()
            sTARGET_TYPE = self.get_argument("TARGET_TYPE").strip()
            sASSAY_TYPE = self.get_argument("ASSAY_TYPE").strip()
            sACTION_TYPE = self.get_argument("ACTION_TYPE").strip()
            sBatches = self.get_argument("batches").strip()
        except:
            logging.error(f"chembl error")
            return
        '''
        if len(sProject) > 2 and len(sELN) > 2:
            sBatches = ''
        '''
        saBatches = sBatches.split()
        sMol = ''
        random_number = str(random.randint(0, 1000000))
        dir_name = f'dist/export/{random_number}'
        os.makedirs(dir_name, exist_ok=True)
        file_path = dir_name + "/COMPOUND_RECORD.tsv"
        molfile_path = dir_name + "/COMPOUND_CTAB.sdf"
        
        with open(file_path, 'w') as compound_record_file, open(molfile_path, 'w') as molfile_file:
            sHeader = ('CIDX', 'RIDX', 'COMPOUND_NAME', 'COMPOUND_KEY')

            compound_record_file.write('\t'.join(map(str, sHeader)) + '\n')
            if len(saBatches) > 10000000:
                chembl_export.exportFromBatches(cur, saBatches, sRIDX, compound_record_file, molfile_file)
            elif len(sProject) > 2 and len(sELN) > 2:
                saCompounds = saBatches
                retVal = chembl_export.exportFromElnProject(cur,
                                                            saCompounds,
                                                            sProject,
                                                            sELN,
                                                            sRIDX,
                                                            sAIDX,
                                                            sTARGET_TYPE,
                                                            sASSAY_TYPE,
                                                            sACTION_TYPE,
                                                            compound_record_file,
                                                            molfile_file,
                                                            dir_name,
                                                            logging)
                if retVal == False:
                    self.write(json.dumps(f'''The combination {sProject} and {sELN} is not valid'''))
                    return

        assay_tsv_name = dir_name + "/ASSAY.tsv"
        with open(assay_tsv_name, 'w') as assay_tsv_file:
            sHeader = (
                'AIDX',
                'RIDX',
                'ASSAY_DESCRIPTION',
                'ASSAY_TYPE',
                'ASSAY_ORGANISM',
                'ASSAY_STRAIN',
                'ASSAY_TAX_ID',
                'ASSAY_TISSUE',
                'ASSAY_CELL_TYPE',
                'ASSAY_SUBCELLULAR_FRACTION',
                'ASSAY_SOURCE',
                'TARGET_TYPE',
                'TARGET_NAME',
                'TARGET_ACCESSION',
                'TARGET_ORGANISM',
                'TARGET_TAX_ID')
            # sAIDX =
            # sRIDX =
            sASSAY_DESCRIPTION = ''
            # sASSAY_TYPE =
            sASSAY_ORGANISM = ''
            sASSAY_STRAIN = ''
            sASSAY_TAX_ID = ''
            sASSAY_TISSUE = ''
            sASSAY_CELL_TYPE = ''
            sASSAY_SUBCELLULAR_FRACTION = ''
            sASSAY_SOURCE = ''
            # sTARGET_TYPE = 
            sTARGET_NAME = ''
            sTARGET_ACCESSION = ''
            sTARGET_ORGANISM = ''
            sTARGET_TAX_ID = ''
            assay_tsv_file.write('\t'.join(map(str, sHeader)) + '\n')
            assay_tsv_file.write(f'''{sAIDX}\t{sRIDX}\t{sASSAY_DESCRIPTION}\t{sASSAY_TYPE}\t{sASSAY_ORGANISM}\t{sASSAY_STRAIN}\t{sASSAY_TAX_ID}\t{sASSAY_TISSUE}\t{sASSAY_CELL_TYPE}\t{sASSAY_SUBCELLULAR_FRACTION}\t{sASSAY_SOURCE}\t{sTARGET_TYPE}\t{sTARGET_NAME}\t{sTARGET_ACCESSION}\t{sTARGET_ORGANISM}\t{sTARGET_TAX_ID}\n''')
        

        reference_tsv_name = dir_name + "/REFERENCE.tsv"
        with open(reference_tsv_name, 'w') as reference_tsv_file:
            sHeader = (
                'RIDX',
                'PUBMED_ID',
                'JOURNAL_NAME',
                'YEAR',
                'VOLUME',
                'ISSUE',
                'FIRST_PAGE',
                'LAST_PAGE',
                'REF_TYPE',
                'TITLE',
                'DOI',
                'PATENT_ID',
                'ABSTRACT',
                'AUTHORS')
            reference_tsv_file.write('\t'.join(map(str, sHeader)) + '\n')
            reference_tsv_file.write(f'''{sRIDX}\t\t\t\t\t\t\t\t\t\t\t\t\t\n''')


        assay_param_tsv_name = dir_name + "/ASSAY_PARAM.tsv"
        with open(assay_param_tsv_name, 'w') as assay_param_tsv_file:
            sHeader = (
                'AIDX',
                'TYPE',
                'RELATION',
                'VALUE',
                'UNITS',
                'TEXT_VALUE',
                'COMMENTS')

            assay_param_tsv_file.write('\t'.join(map(str, sHeader)) + '\n')
            assay_param_tsv_file.write(f'''{sAIDX}\t\t\t\t\t\t\n''')

            
        zip_filepath = os.path.join(dir_name, "ChEMBL_deposit.zip")

        with zipfile.ZipFile(zip_filepath, 'w') as zipf:
            zipf.write(os.path.join(dir_name, "COMPOUND_RECORD.tsv"), arcname="COMPOUND_RECORD.tsv")
            zipf.write(os.path.join(dir_name, "COMPOUND_CTAB.sdf"), arcname="COMPOUND_CTAB.sdf")
            zipf.write(os.path.join(dir_name, "ACTIVITY.tsv"), arcname="ACTIVITY.tsv")
            zipf.write(os.path.join(dir_name, "ASSAY.tsv"), arcname="ASSAY.tsv")
            zipf.write(os.path.join(dir_name, "REFERENCE.tsv"), arcname="REFERENCE.tsv")
            zipf.write(os.path.join(dir_name, "ASSAY_PARAM.tsv"), arcname="ASSAY_PARAM.tsv")

        
        self.write(json.dumps(f'''<a href=https://esox3.scilifelab.se/chemreg/dist/export/{random_number}/ChEMBL_deposit.zip>ChEMBL_deposit.zip</a>'''))
