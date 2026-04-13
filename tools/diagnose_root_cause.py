#!/usr/bin/env python3
"""
Root cause diagnostic: trace where wrong smiles_std came from.
Uses the same DB connection as the server (via config.py).
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import MySQLdb
import config
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

print(f"RDKit version: {Chem.rdBase.rdkitVersion}")

conn = MySQLdb.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password'],
    db='bcpvs',
    charset='utf8mb4',
    use_unicode=True
)
conn.autocommit(True)
cur = conn.cursor()


def rdkit_standardize_molfile(molfile):
    """Reproduce cleanStructureRDKit's pipeline (RDKit part only)."""
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return None
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    try:
        clean_mol = rdMolStandardize.Cleanup(mol, params)
        parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
        uncharged = rdMolStandardize.Uncharger().uncharge(parent_mol)
        te = rdMolStandardize.TautomerEnumerator(params)
        canon = te.Canonicalize(uncharged)
        return Chem.MolToSmiles(canon)
    except Exception as e:
        return f"ERROR: {e}"


# The 11 affected compounds from the SDF double-registration
AFFECTED = [
    'CBK028998', 'CBK228517', 'CBK278468', 'CBK278635',
    'CBK291635', 'CBK291792', 'CBK389719', 'CBK699880',
    'CBK714175', 'CBK714200',
]

print("=" * 110)
print("PART 1: For each affected compound, compare stored smiles_std vs RDKit from chem_reg MOLFILE vs JCMOL")
print("=" * 110)

for cid in AFFECTED:
    print(f"\n{'='*80}")
    print(f"COMPOUND: {cid}")
    print(f"{'='*80}")

    # 1. Get stored smiles_std
    cur.execute("SELECT smiles_std, created_date FROM bcpvs.compound WHERE compound_id = %s", (cid,))
    row = cur.fetchone()
    if not row:
        print("  NOT FOUND in compound table")
        continue
    stored_smiles, created_date = row
    print(f"  compound.smiles_std  = {stored_smiles}")
    print(f"  created_date         = {created_date}")

    # 2. Get ALL batches for this compound
    cur.execute("""SELECT notebook_ref, chemreg_regno, chemspec_regno, submittal_date
                   FROM bcpvs.batch WHERE compound_id = %s ORDER BY submittal_date""", (cid,))
    batches = cur.fetchall()
    print(f"  Batches: {len(batches)}")
    for (nb, cr, cs, sd) in batches:
        print(f"    {nb} chemreg_regno={cr} chemspec_regno={cs} date={sd}")

    # 3. Get earliest chem_reg molfile
    cur.execute("""SELECT ci.MOLFILE, ci.regno, ci.RDATE, ci.JPAGE
                   FROM chem_reg.chem_info ci
                   WHERE ci.COMPOUND_ID = %s AND ci.MOLFILE IS NOT NULL
                   ORDER BY ci.regno ASC""", (cid,))
    ci_rows = cur.fetchall()
    for (molfile, regno, rdate, jpage) in ci_rows:
        rdkit_smi = rdkit_standardize_molfile(molfile) if molfile else None
        match_stored = (rdkit_smi == stored_smiles)
        print(f"    chem_reg: regno={regno} jpage={jpage} rdate={rdate}")
        print(f"      RDKit from MOLFILE = {rdkit_smi}")
        print(f"      Matches stored?    = {match_stored}")

    # 4. Get JCMOL mol
    cur.execute("SELECT mol FROM bcpvs.JCMOL_MOLTABLE WHERE compound_id = %s", (cid,))
    jrow = cur.fetchone()
    if jrow and jrow[0]:
        jcmol = jrow[0]
        rdkit_jcmol = rdkit_standardize_molfile(jcmol)
        print(f"  RDKit from JCMOL     = {rdkit_jcmol}")
        print(f"  JCMOL matches stored?= {rdkit_jcmol == stored_smiles}")
    else:
        print(f"  JCMOL: NOT FOUND")

    # 5. Molcart round-trip: send stored SMILES through mol2bin/bin2mol
    try:
        cur.execute("SELECT CONVERT(bin2mol(mol2bin(%s, 'smiles')) USING UTF8)", (stored_smiles,))
        mc_row = cur.fetchone()
        if mc_row and mc_row[0]:
            molcart_molfile = mc_row[0]
            rdkit_from_molcart = rdkit_standardize_molfile(molcart_molfile)
            print(f"  Molcart roundtrip    = {rdkit_from_molcart}")
            print(f"  Roundtrip matches?   = {rdkit_from_molcart == stored_smiles}")
        else:
            print(f"  Molcart roundtrip: NULL result")
    except Exception as e:
        print(f"  Molcart roundtrip ERROR: {e}")


print("\n\n")
print("=" * 110)
print("PART 2: Check smiles_std column collation (case-sensitivity)")
print("=" * 110)

cur.execute("""SELECT COLUMN_NAME, CHARACTER_SET_NAME, COLLATION_NAME
               FROM INFORMATION_SCHEMA.COLUMNS
               WHERE TABLE_SCHEMA = 'bcpvs'
               AND TABLE_NAME = 'compound'
               AND COLUMN_NAME IN ('smiles_std', 'smiles_std_string')""")
for (col, charset, coll) in cur.fetchall():
    is_ci = '_ci' in (coll or '')
    print(f"  {col}: charset={charset}, collation={coll}  {'*** CASE-INSENSITIVE ***' if is_ci else 'case-sensitive'}")

# Direct test: do two compounds exist that differ only in case?
cur.execute("""SELECT compound_id, smiles_std FROM bcpvs.compound WHERE compound_id = 'CBK278468'""")
row = cur.fetchone()
if row:
    wrong = row[1]  # NC(=O)c1ccncc1 (pyridine - aromatic)
    # Search with the saturated version
    saturated = 'NC(=O)C1CCNCC1'
    cur.execute("SELECT compound_id FROM bcpvs.compound WHERE smiles_std = %s", (saturated,))
    found = cur.fetchall()
    print(f"\n  CBK278468 stored: '{wrong}'")
    print(f"  Search for '{saturated}': found {[r[0] for r in found]}")
    cur.execute("SELECT compound_id FROM bcpvs.compound WHERE smiles_std = %s", (wrong,))
    found2 = cur.fetchall()
    print(f"  Search for '{wrong}': found {[r[0] for r in found2]}")


print("\n\n")
print("=" * 110)
print("PART 3: Check the duplicate compound_ids (from 2nd SDF file registration)")
print("=" * 110)

# These are compounds that were created as duplicates when registering the second SDF file
# They should have the CORRECT SMILES (since RDKit produces correct results)
# Find batches with notebook_ref matching the EH833xxxx pattern from the CSV
cur.execute("""SELECT b.compound_id, b.notebook_ref, b.chemreg_regno, c.smiles_std, c.created_date
               FROM bcpvs.batch b
               JOIN bcpvs.compound c ON b.compound_id = c.compound_id
               WHERE b.notebook_ref LIKE 'EH833%'
               AND b.compound_id LIKE 'CBK7141%'
               ORDER BY b.compound_id""")
for row in cur.fetchall():
    print(f"  {row}")

# Let's look at what compound_ids were created around the same time as the affected compounds
print("\n  Recently created compounds (should show both originals and duplicates):")
for cid in AFFECTED:
    cur.execute("""SELECT c.compound_id, c.smiles_std, c.created_date
                   FROM bcpvs.compound c 
                   WHERE c.smiles_std IN (
                       SELECT c2.smiles_std FROM bcpvs.compound c2 WHERE c2.compound_id = %s
                   )""", (cid,))
    rows = cur.fetchall()
    if len(rows) > 1:
        print(f"\n  Multiple compounds with same smiles_std as {cid}:")
        for r in rows:
            print(f"    {r}")


print("\n\n")
print("=" * 110)
print("PART 4: Molcart behavior test - does molcart change aromaticity?")
print("=" * 110)

test_smiles = [
    ('NC(=O)C1CCNCC1', 'Piperidine amide (saturated)'),
    ('NC(=O)c1ccncc1', 'Nicotinamide (aromatic)'),
    ('O=C(NC1CCCCC1)c1ccccc1', 'N-cyclohexylbenzamide'),
    ('O=C(Nc1ccccc1)c1ccccc1', 'Benzanilide'),
    ('CNC(=O)c1ccncc1', 'N-methyl nicotinamide (aromatic)'),
    ('CNC(=O)C1CCCNC1', 'N-methyl piperidine-3-carboxamide (saturated)'),
]

for smi, name in test_smiles:
    try:
        # mol2bin → bin2mol roundtrip (same as cleanStructureRDKit)
        cur.execute("SELECT CONVERT(bin2mol(mol2bin(%s, 'smiles')) USING UTF8)", (smi,))
        mc_row = cur.fetchone()
        if mc_row and mc_row[0]:
            mc_molfile = mc_row[0]
            # Process the molcart molfile through RDKit
            mc_rdkit = rdkit_standardize_molfile(mc_molfile)
            match = (mc_rdkit == smi)
            print(f"  {name}")
            print(f"    Input:     {smi}")
            print(f"    Molcart→RDKit: {mc_rdkit}")
            print(f"    Match: {match}")
            if not match:
                print(f"    *** MOLCART CHANGED THE STRUCTURE! ***")
        else:
            print(f"  {name}: molcart returned NULL for {smi}")
    except Exception as e:
        print(f"  {name}: ERROR: {e}")


print("\n\n")
print("=" * 110)
print("PART 5: Full cleanStructureRDKit simulation including molcart roundtrip")
print("=" * 110)

for smi, name in test_smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue
    
    # Standardize with RDKit (same as cleanStructureRDKit)
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    clean_mol = rdMolStandardize.Cleanup(mol, params)
    parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    uncharged = rdMolStandardize.Uncharger().uncharge(parent_mol)
    te = rdMolStandardize.TautomerEnumerator(params)
    canon = te.Canonicalize(uncharged)
    mV4 = Chem.MolToSmiles(canon)
    mV4_escaped = mV4.replace('\\', '\\\\')
    
    # Molcart roundtrip
    try:
        cur.execute(f"SELECT CONVERT(bin2mol(mol2bin('{mV4_escaped}', 'smiles')) USING UTF8)")
        mc_row = cur.fetchone()
        molcart_molfile = mc_row[0] if mc_row else None
    except Exception as e:
        molcart_molfile = None
        print(f"  Molcart error: {e}")
    
    if molcart_molfile:
        # Now process the molcart molfile through RDKit again (this is what gets stored in JCMOL)
        rdkit_from_mc = rdkit_standardize_molfile(molcart_molfile)
        print(f"  {name}:")
        print(f"    Original SMILES:    {smi}")
        print(f"    RDKit stdSMILES:    {mV4}")
        print(f"    Molcart→RDKit:      {rdkit_from_mc}")
        if mV4 != rdkit_from_mc:
            print(f"    *** DIVERGENCE: molcart molfile gives different RDKit SMILES! ***")

cur.close()
conn.close()
