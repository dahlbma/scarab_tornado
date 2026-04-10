"""
Test to run ON THE SERVER (RDKit 2022.3) to determine:
1. What SMILES does cleanStructureRDKit produce from the actual SDF molfiles?
2. Is the tautomer canonicalization non-deterministic across multiple calls?

This tests the actual registration path: molfile → RDKit → tautomer → SMILES
"""
import re
import io
import sys
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

print(f"RDKit version: {Chem.rdBase.rdkitVersion}")
print()

AFFECTED_SUPPLIER_IDS = [
    'Z104477230', 'Z199507048', 'Z27749575', 'Z285662830',
    'Z31407959', 'Z31484244', 'Z31484247', 'Z31484290',
    'Z32016358', 'Z44585993',
]

# Known values from the database
DB_REGISTERED_SMILES = {
    'Z104477230': 'NC(=O)c1ccncc1',       # piperidine→pyridine (wrong)
    'Z199507048': 'O=C(NC1CCNCC1)C1CCCCC1', # benzene→cyclohexane (wrong)
    'Z27749575':  'O=C(NCc1ccccc1)c1ccccc1', # cyclohexane→benzene (wrong)
    'Z285662830': 'CNC(=O)c1cccnc1',       # piperidine→pyridine (wrong)
    'Z31407959':  'O=C(NCC1CCCO1)c1cccc(Cl)c1', # furan→THF (wrong)
    'Z31484244':  'O=C(Nc1ccccc1)c1cccc(O)c1',  # cyclohexane→benzene (wrong)
    'Z31484247':  'O=C(Nc1ccccc1)c1cccnc1',     # cyclohexane→benzene (wrong)
    'Z31484290':  'O=C(Nc1ccccc1)c1ccccc1',     # cyclohexane→benzene (wrong)
    'Z32016358':  'CNC(=O)C1CCNCC1',       # pyridine→piperidine (wrong)
    'Z44585993':  'O=C(NCC1CCCO1)NC1CCCCC1', # benzene→cyclohexane (wrong)
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
    """
    Same pipeline as the server's cleanStructureRDKit, but without the MySQL
    mol2bin/bin2mol step. Returns (molblock_str, stdSMILES).
    """
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
    canonical_tautomer = te.Canonicalize(uncharged)

    stdSMILES = Chem.MolToSmiles(canonical_tautomer)
    return Chem.MolToMolBlock(canonical_tautomer), stdSMILES


def enumerate_tautomers_from_molfile(molfile):
    """Enumerate all tautomers from a molfile."""
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return [], None

    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False

    clean_mol = rdMolStandardize.Cleanup(mol, params)
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(parent_clean_mol)

    # Also get the SMILES before tautomer step
    pre_tautomer_smiles = Chem.MolToSmiles(uncharged)

    te = rdMolStandardize.TautomerEnumerator(params)
    tautomers = te.Enumerate(uncharged)

    return [Chem.MolToSmiles(t) for t in tautomers], pre_tautomer_smiles


# Determine sdf_dir based on where the script is run
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
# Try to find the SDF files
for candidate in [os.path.join(parent_dir, 'chemReg'), '.', '..', '../chemReg']:
    test_file = os.path.join(candidate, 'CBCS0523_CC04606_208cpds_Carlsson.sdf')
    if os.path.exists(test_file):
        sdf_dir = candidate
        break
else:
    print("ERROR: Cannot find SDF files. Please run from the tools/ or chemReg/ directory.")
    sys.exit(1)

print(f"Reading SDF files from: {os.path.abspath(sdf_dir)}")
file1 = os.path.join(sdf_dir, 'CBCS0523_CC04606_208cpds_Carlsson.sdf')
file2 = os.path.join(sdf_dir, 'CBCS0523_CC04640_192cpds_Carlsson.sdf')

mols1 = parse_sdf(file1)
mols2 = parse_sdf(file2)

print(f"File 1: {len(mols1)} molecules")
print(f"File 2: {len(mols2)} molecules")

print()
print("=" * 100)
print("TEST 1: What SMILES does cleanStructureRDKit produce from the molfile?")
print("        Compare with what was stored in the database.")
print("=" * 100)

for sid in AFFECTED_SUPPLIER_IDS:
    molblock = mols1.get(sid)
    if not molblock:
        print(f"\n  {sid}: NOT FOUND in File 1")
        continue

    _, std_smi = cleanStructureRDKit_local(molblock)
    db_smi = DB_REGISTERED_SMILES.get(sid, '???')

    # Get raw SMILES (no standardization)
    mol = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=True)
    raw_smi = Chem.MolToSmiles(mol) if mol else "PARSE_FAILED"

    print(f"\n--- {sid} ---")
    print(f"  Raw SMILES (from molfile, no standardization): {raw_smi}")
    print(f"  Standardized SMILES (this run):                {std_smi}")
    print(f"  Database registered SMILES:                    {db_smi}")

    if std_smi == db_smi:
        print(f"  → Current result MATCHES database (same wrong form)")
    elif std_smi == raw_smi:
        print(f"  → Current result matches raw = CORRECT (tautomer did not corrupt)")
        print(f"  → But database has DIFFERENT wrong form!")
        print(f"  → THIS PROVES NON-DETERMINISM between runs")
    else:
        print(f"  → Current result differs from BOTH raw and database!")


print()
print("=" * 100)
print("TEST 2: Tautomer enumeration from molfiles")
print("        Check if aromatic/saturated forms appear as tautomers")
print("=" * 100)

for sid in AFFECTED_SUPPLIER_IDS:
    molblock = mols1.get(sid)
    if not molblock:
        continue

    tautomers, pre_taut_smi = enumerate_tautomers_from_molfile(molblock)

    print(f"\n--- {sid} ---")
    print(f"  Pre-tautomer SMILES: {pre_taut_smi}")
    print(f"  Number of tautomers: {len(tautomers)}")
    for i, t in enumerate(tautomers):
        print(f"    Tautomer {i}: {t}")

    # Check if any tautomer has different aromaticity
    mol_pre = Chem.MolFromSmiles(pre_taut_smi)
    pre_arom = sum(1 for a in mol_pre.GetAtoms() if a.GetIsAromatic()) if mol_pre else -1

    for i, t in enumerate(tautomers):
        mol_t = Chem.MolFromSmiles(t)
        t_arom = sum(1 for a in mol_t.GetAtoms() if a.GetIsAromatic()) if mol_t else -1
        if t_arom != pre_arom:
            print(f"    *** Tautomer {i} has DIFFERENT aromaticity! ({pre_arom} → {t_arom} aromatic atoms)")
            print(f"    *** This is the RDKit bug: '{pre_taut_smi}' <-> '{t}'")


print()
print("=" * 100)
print("TEST 3: Non-determinism check - run standardization 50 times")
print("=" * 100)

for sid in AFFECTED_SUPPLIER_IDS:
    molblock = mols1.get(sid)
    if not molblock:
        continue

    results = set()
    for _ in range(50):
        _, smi = cleanStructureRDKit_local(molblock)
        if smi:
            results.add(smi)

    if len(results) > 1:
        print(f"  {sid}: NON-DETERMINISTIC! Got {len(results)} different SMILES:")
        for r in sorted(results):
            print(f"    - {r}")
    else:
        smi = list(results)[0] if results else None
        db_smi = DB_REGISTERED_SMILES.get(sid, '???')
        if smi == db_smi:
            match_str = "(matches DB)"
        else:
            match_str = f"(DB has: {db_smi}) - DIFFERS FROM DB!"
        print(f"  {sid}: Deterministic within this run: {smi} {match_str}")

print()
print("NOTE: If all are 'Deterministic within this run' but some 'DIFFER FROM DB',")
print("      this means the non-determinism occurs ACROSS process restarts.")
print("      Try restarting and running this script again to see if results change.")
