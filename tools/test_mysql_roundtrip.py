"""
Test to run ON THE SERVER to determine if the MySQL round-trip of the molfile
causes cleanStructureRDKit to produce different (wrong) SMILES.

HYPOTHESIS: The molfile stored via f-string SQL interpolation in ChemRegAddMol:
    sSql = f"insert into ... values (..., '{molfile}')"
gets modified by MySQL (e.g. backslash escapes), and when BcpvsRegCompound
reads it back, cleanStructureRDKit produces a different SMILES.

This test:
1. Reads original molblocks from the SDF files (as getNextMolecule would)
2. Simulates the MySQL round-trip by inserting/reading from a temp table
3. Compares cleanStructureRDKit output from original vs round-tripped molfile
4. Also reads the ACTUAL stored molfile from chem_info for the affected regnos

REQUIRES: Run from the scarab_tornado directory (needs config.py and mydb.py)
"""
import re
import io
import sys
import os

# Add parent dir to path so we can import config and mydb
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, parent_dir)

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import mydb

print(f"RDKit version: {Chem.rdBase.rdkitVersion}")
print()

# Connect to the database using the same connection as the server
db = mydb.disconnectSafeConnect()
cur = db.cursor()

AFFECTED = {
    'Z104477230': {'db_smiles': 'NC(=O)c1ccncc1'},
    'Z199507048': {'db_smiles': 'O=C(NC1CCNCC1)C1CCCCC1'},
    'Z27749575':  {'db_smiles': 'O=C(NCc1ccccc1)c1ccccc1'},
    'Z285662830': {'db_smiles': 'CNC(=O)c1cccnc1'},
    'Z31407959':  {'db_smiles': 'O=C(NCC1CCCO1)c1cccc(Cl)c1'},
    'Z31484244':  {'db_smiles': 'O=C(Nc1ccccc1)c1cccc(O)c1'},
    'Z31484247':  {'db_smiles': 'O=C(Nc1ccccc1)c1cccnc1'},
    'Z31484290':  {'db_smiles': 'O=C(Nc1ccccc1)c1ccccc1'},
    'Z32016358':  {'db_smiles': 'CNC(=O)C1CCNCC1'},
    'Z44585993':  {'db_smiles': 'O=C(NCC1CCCO1)NC1CCCCC1'},
}


def parse_sdf(filename):
    """Simulate sdfreg.py's getNextMolecule exactly."""
    results = {}
    with open(filename, 'rb') as f:
        content = f.read()
    f = io.BytesIO(content)
    while True:
        sMol = b""
        iCount = 0
        found_end = False
        while True:
            iCount += 1
            if iCount > 10000:
                break
            line = f.readline()
            if not line:
                break
            line = line.replace(b'\r\n', b'\n')
            line = line.replace(b"'", b"")
            if iCount == 1:
                line = b'id' + line
            sMol += line
            if b'$$$$' in line:
                found_end = True
                break
        if not found_end or len(sMol) < 4:
            break
        mol_str = sMol.decode(errors='replace')
        pattern = r'>\s*<External ID>\s*\n(\S+)'
        match = re.search(pattern, mol_str)
        if match:
            ext_id = match.group(1).strip()
        else:
            continue
        m_end_idx = mol_str.find('M  END')
        if m_end_idx >= 0:
            molblock = mol_str[:m_end_idx + 6]
        else:
            continue
        results[ext_id] = molblock
    return results


def cleanStructureRDKit_local(molfile):
    """Same as server, without the mol2bin/bin2mol step."""
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return None, None
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    try:
        clean_mol = rdMolStandardize.Cleanup(mol, params)
    except:
        return None, None
    try:
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    except:
        return None, None
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(parent_clean_mol)
    te = rdMolStandardize.TautomerEnumerator(params)
    canonical = te.Canonicalize(uncharged)
    stdSMILES = Chem.MolToSmiles(canonical)
    return Chem.MolToMolBlock(canonical), stdSMILES


def cleanStructureRDKit_with_molcart(molfile):
    """Same as server's cleanStructureRDKit, INCLUDING the mol2bin/bin2mol step."""
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return None, None
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    try:
        clean_mol = rdMolStandardize.Cleanup(mol, params)
    except:
        return None, None
    try:
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    except:
        return None, None
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(parent_clean_mol)
    te = rdMolStandardize.TautomerEnumerator(params)
    canonical = te.Canonicalize(uncharged)
    mV4 = Chem.MolToSmiles(canonical)
    mV4 = mV4.replace('\\', '\\\\')

    # This is the Molcart round-trip that the server does
    sSql = f'''select CONVERT(bin2mol(mol2bin('{mV4}' , 'smiles')) USING UTF8)'''
    cur.execute(sSql)
    saRes = cur.fetchall()
    if len(saRes) == 1:
        return saRes[0][0], mV4
    else:
        return None, None


# Find SDF files
for candidate in [os.path.join(parent_dir, 'chemReg'), '.', '..', '../chemReg']:
    test_file = os.path.join(candidate, 'CBCS0523_CC04606_208cpds_Carlsson.sdf')
    if os.path.exists(test_file):
        sdf_dir = candidate
        break
else:
    print("ERROR: Cannot find SDF files.")
    sys.exit(1)

print(f"Reading SDF files from: {os.path.abspath(sdf_dir)}")
file1 = os.path.join(sdf_dir, 'CBCS0523_CC04606_208cpds_Carlsson.sdf')
mols1 = parse_sdf(file1)
print(f"File 1: {len(mols1)} molecules")

# ============================================================================
print()
print("=" * 100)
print("TEST A: MySQL molfile round-trip")
print("         Insert original molfile into temp table, read it back, compare")
print("=" * 100)

# Create a temporary table
try:
    cur.execute("DROP TABLE IF EXISTS chem_reg.temp_molfile_test")
    cur.execute("""CREATE TABLE chem_reg.temp_molfile_test (
        id INT PRIMARY KEY AUTO_INCREMENT,
        ext_id VARCHAR(50),
        molfile MEDIUMTEXT
    )""")
except Exception as e:
    print(f"Could not create temp table in chem_reg: {e}")
    print("Trying chem_reg_test...")
    try:
        cur.execute("DROP TABLE IF EXISTS chem_reg_test.temp_molfile_test")
        cur.execute("""CREATE TABLE chem_reg_test.temp_molfile_test (
            id INT PRIMARY KEY AUTO_INCREMENT,
            ext_id VARCHAR(50),
            molfile MEDIUMTEXT
        )""")
    except Exception as e2:
        print(f"Could not create temp table: {e2}")
        print("SKIPPING MySQL round-trip test")
        sys.exit(1)

test_db = 'chem_reg'  # or 'chem_reg_test'

for sid in sorted(AFFECTED.keys()):
    molblock = mols1.get(sid)
    if not molblock:
        continue

    # The full mol content including tags, as sent by sdfreg.py client
    # The molfile would be decoded as latin-1 in the client:
    # dTags['molfile'] = sMol.decode('latin-1')
    # Then sent via HTTP POST and stored via f-string SQL

    # Insert using f-string SQL interpolation (same as the server does)
    try:
        sSql = f"""INSERT INTO {test_db}.temp_molfile_test (ext_id, molfile)
                   VALUES ('{sid}', '{molblock}')"""
        cur.execute(sSql)
    except Exception as e:
        print(f"  {sid}: INSERT FAILED: {e}")
        continue

    # Read it back
    sSql = f"SELECT molfile FROM {test_db}.temp_molfile_test WHERE ext_id = '{sid}'"
    cur.execute(sSql)
    row = cur.fetchall()
    if not row:
        print(f"  {sid}: READ BACK FAILED")
        continue
    molblock_from_db = row[0][0]

    # Compare
    if molblock == molblock_from_db:
        print(f"  {sid}: Molfile IDENTICAL after MySQL round-trip")
    else:
        print(f"  {sid}: *** Molfile CHANGED after MySQL round-trip! ***")
        # Find differences
        lines_orig = molblock.split('\n')
        lines_db = molblock_from_db.split('\n')
        print(f"    Original: {len(lines_orig)} lines, DB: {len(lines_db)} lines")
        for i, (l1, l2) in enumerate(zip(lines_orig, lines_db)):
            if l1 != l2:
                print(f"    Line {i} differs:")
                print(f"      Original: {repr(l1)}")
                print(f"      From DB:  {repr(l2)}")
                if i > 5:
                    print(f"    ... (showing first 5 differences)")
                    break

        # Process both through cleanStructureRDKit
        _, smi_orig = cleanStructureRDKit_local(molblock)
        _, smi_db = cleanStructureRDKit_local(molblock_from_db)
        print(f"    SMILES from original: {smi_orig}")
        print(f"    SMILES from DB copy:  {smi_db}")
        if smi_orig != smi_db:
            print(f"    *** SMILES DIFFER! This explains the bug! ***")

# Cleanup temp table
try:
    cur.execute(f"DROP TABLE IF EXISTS {test_db}.temp_molfile_test")
except:
    pass

# ============================================================================
print()
print("=" * 100)
print("TEST B: Read ACTUAL stored molfiles from chem_info for affected regnos")
print("         Compare with SDF originals")
print("=" * 100)

# Find regnos for the affected external_ids
for sid in sorted(AFFECTED.keys()):
    molblock_orig = mols1.get(sid)
    if not molblock_orig:
        continue

    # Look up in chem_info by external_id
    sSql = f"""SELECT regno, molfile, compound_id, jpage
               FROM chem_reg.chem_info
               WHERE external_id = '{sid}'
               ORDER BY regno"""
    cur.execute(sSql)
    rows = cur.fetchall()

    if not rows:
        # Try test database
        sSql = f"""SELECT regno, molfile, compound_id, jpage
                   FROM chem_reg_test.chem_info
                   WHERE external_id = '{sid}'
                   ORDER BY regno"""
        cur.execute(sSql)
        rows = cur.fetchall()

    if not rows:
        print(f"\n  {sid}: Not found in chem_info")
        continue

    print(f"\n--- {sid} ({len(rows)} registration(s)) ---")
    for (regno, molfile_db, compound_id, jpage) in rows:
        print(f"  Regno: {regno}, Compound: {compound_id}, JPage: {jpage}")

        if molfile_db is None:
            print(f"    Molfile is NULL in database!")
            continue

        # Extract just the molblock (up to M  END) from the DB molfile
        m_end_idx = molfile_db.find('M  END')
        if m_end_idx >= 0:
            molblock_db = molfile_db[:m_end_idx + 6]
        else:
            molblock_db = molfile_db

        # Compare molblocks
        if molblock_orig == molblock_db:
            print(f"    Molblock: IDENTICAL to SDF original")
        else:
            print(f"    Molblock: *** DIFFERS from SDF original! ***")
            lines_orig = molblock_orig.split('\n')
            lines_db = molblock_db.split('\n')
            print(f"      Original: {len(lines_orig)} lines, DB: {len(lines_db)} lines")
            diff_count = 0
            for i, (l1, l2) in enumerate(zip(lines_orig, lines_db)):
                if l1 != l2:
                    diff_count += 1
                    if diff_count <= 5:
                        print(f"      Line {i}: orig={repr(l1)} db={repr(l2)}")
            if diff_count > 5:
                print(f"      ... ({diff_count} total differing lines)")
            if len(lines_orig) != len(lines_db):
                print(f"      Line count differs: {len(lines_orig)} vs {len(lines_db)}")

        # Process both through cleanStructureRDKit and compare
        _, smi_orig = cleanStructureRDKit_local(molblock_orig)
        _, smi_db = cleanStructureRDKit_local(molblock_db)
        db_smiles = AFFECTED[sid]['db_smiles']

        print(f"    SMILES from SDF original:  {smi_orig}")
        print(f"    SMILES from DB molfile:    {smi_db}")
        print(f"    SMILES in compound table:  {db_smiles}")

        if smi_orig != smi_db:
            print(f"    *** DIFFERENT SMILES from original vs DB! MySQL modified the molfile! ***")
        elif smi_db != db_smiles:
            print(f"    *** DB molfile gives correct SMILES, but compound table has wrong! ***")
            print(f"    → The wrong SMILES was produced during BcpvsRegCompound under different conditions")


# ============================================================================
print()
print("=" * 100)
print("TEST C: Check if Molcart mol2bin/bin2mol changes the molecule")
print("         (the last step of cleanStructureRDKit)")
print("=" * 100)

for sid in sorted(AFFECTED.keys()):
    molblock = mols1.get(sid)
    if not molblock:
        continue

    # Full pipeline with Molcart
    molfile_molcart, smi_with_molcart = cleanStructureRDKit_with_molcart(molblock)
    # Pipeline without Molcart
    _, smi_without_molcart = cleanStructureRDKit_local(molblock)

    print(f"  {sid}:")
    print(f"    SMILES (RDKit only):          {smi_without_molcart}")
    print(f"    SMILES (RDKit+Molcart):       {smi_with_molcart}")

    if smi_with_molcart != smi_without_molcart:
        print(f"    *** Molcart changed the SMILES! ***")
    else:
        print(f"    Molcart did not change SMILES")

    # Also check: does the Molcart-returned molfile produce different SMILES?
    if molfile_molcart:
        mol_from_molcart = Chem.MolFromMolBlock(molfile_molcart, removeHs=False, sanitize=True)
        if mol_from_molcart:
            smi_from_molcart_molfile = Chem.MolToSmiles(mol_from_molcart)
            if smi_from_molcart_molfile != smi_without_molcart:
                print(f"    Molcart molfile SMILES:        {smi_from_molcart_molfile}")
                print(f"    *** Molcart molfile gives DIFFERENT molecule! ***")

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)
print("""
If TEST A shows molfile changed by MySQL → MySQL backslash escape is the root cause
If TEST B shows DB molfile differs from original → confirms MySQL modification
If TEST C shows Molcart changes SMILES → Molcart is involved in the bug

If NONE of the above show differences → non-determinism occurred at registration
time under different process conditions (PYTHONHASHSEED, server version, etc.)
""")
