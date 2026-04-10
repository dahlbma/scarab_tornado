"""
Final diagnostic test to identify why the wrong SMILES ended up in the compound table
and why duplicate detection failed.

WHAT WE KNOW:
- Current RDKit 2022.3 produces CORRECT SMILES from the molfiles
- The compound table has WRONG smiles_std (aromatic/saturated swaps)
- The wrong SMILES must have been produced by a DIFFERENT RDKit/environment
- RDKit 2022.3 is deterministic

THIS TEST CHECKS:
1. MySQL column collation - is smiles_std comparison case-sensitive?
   (In SMILES, c=aromatic, C=aliphatic - case matters!)
2. When were these compounds first registered? Was it by the current server?
3. Are there other RDKit installations that might have been used?
4. Exactly which step in BcpvsRegCompound produced the wrong SMILES?

REQUIRES: Run from the scarab_tornado directory (needs config.py and mydb.py)
"""
import sys
import os
import subprocess

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, parent_dir)

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import mydb

print(f"RDKit version: {Chem.rdBase.rdkitVersion}")
print(f"Python: {sys.executable}")
print()

db = mydb.disconnectSafeConnect()
cur = db.cursor()

# The affected compound_ids and their known wrong smiles_std
AFFECTED_COMPOUNDS = {
    'CBK278468': {'wrong_smi': 'NC(=O)c1ccncc1',       'correct_smi': 'NC(=O)C1CCNCC1'},
    'CBK389719': {'wrong_smi': 'O=C(NC1CCNCC1)C1CCCCC1', 'correct_smi': 'O=C(NC1CCNCC1)c1ccccc1'},
    'CBK699880': {'wrong_smi': 'O=C(NCc1ccccc1)c1ccccc1', 'correct_smi': 'O=C(NCc1ccccc1)C1CCCCC1'},
    'CBK714175': {'wrong_smi': 'CNC(=O)c1cccnc1',       'correct_smi': 'CNC(=O)C1CCCNC1'},
    'CBK714200': {'wrong_smi': 'O=C(NCC1CCCO1)c1cccc(Cl)c1', 'correct_smi': 'O=C(NCc1ccco1)c1cccc(Cl)c1'},
    'CBK278635': {'wrong_smi': 'O=C(Nc1ccccc1)c1cccc(O)c1', 'correct_smi': 'O=C(NC1CCCCC1)c1cccc(O)c1'},
    'CBK291635': {'wrong_smi': 'O=C(Nc1ccccc1)c1cccnc1',   'correct_smi': 'O=C(NC1CCCCC1)c1cccnc1'},
    'CBK228517': {'wrong_smi': 'O=C(Nc1ccccc1)c1ccccc1',   'correct_smi': 'O=C(NC1CCCCC1)c1ccccc1'},
    'CBK291792': {'wrong_smi': 'CNC(=O)C1CCNCC1',       'correct_smi': 'CNC(=O)c1ccncc1'},
    'CBK028998': {'wrong_smi': 'O=C(NCC1CCCO1)NC1CCCCC1', 'correct_smi': 'O=C(NCC1CCCO1)Nc1ccccc1'},
}


print("=" * 100)
print("TEST 1: MySQL collation check for smiles_std column")
print("        Is string comparison case-sensitive for SMILES?")
print("=" * 100)

for db_name in ['bcpvs', 'bcpvs_test']:
    try:
        sSql = f"""SELECT COLUMN_NAME, CHARACTER_SET_NAME, COLLATION_NAME
                   FROM INFORMATION_SCHEMA.COLUMNS
                   WHERE TABLE_SCHEMA = '{db_name}'
                   AND TABLE_NAME = 'compound'
                   AND COLUMN_NAME IN ('smiles_std', 'smiles_std_string', 'compound_id')"""
        cur.execute(sSql)
        rows = cur.fetchall()
        if rows:
            print(f"\n  Database: {db_name}")
            for (col, charset, collation) in rows:
                is_ci = '_ci' in (collation or '')
                print(f"    {col}: charset={charset}, collation={collation}")
                if is_ci:
                    print(f"      *** CASE-INSENSITIVE! ***")
                    print(f"      This means 'c1ccccc1' WOULD match 'C1CCCCC1' in WHERE clause")
                else:
                    print(f"      Case-sensitive comparison")
    except Exception as e:
        print(f"\n  Database {db_name}: {e}")

# Direct test: does a case-different SMILES match?
print("\n  Direct test with actual data:")
for db_name in ['bcpvs', 'bcpvs_test']:
    try:
        cid = 'CBK278468'
        wrong_smi = AFFECTED_COMPOUNDS[cid]['wrong_smi']
        correct_smi = AFFECTED_COMPOUNDS[cid]['correct_smi']

        sSql = f"SELECT compound_id, smiles_std FROM {db_name}.compound WHERE compound_id = '{cid}'"
        cur.execute(sSql)
        rows = cur.fetchall()
        if rows:
            stored_smi = rows[0][1]
            print(f"\n  {db_name}.compound: {cid} has smiles_std = '{stored_smi}'")

            # Try searching with the correct SMILES
            sSql = f"SELECT compound_id FROM {db_name}.compound WHERE smiles_std = '{correct_smi}'"
            cur.execute(sSql)
            found = cur.fetchall()
            print(f"  Search with correct '{correct_smi}': found {len(found)} results: {[r[0] for r in found]}")

            # Try searching with the wrong SMILES
            sSql = f"SELECT compound_id FROM {db_name}.compound WHERE smiles_std = '{wrong_smi}'"
            cur.execute(sSql)
            found = cur.fetchall()
            print(f"  Search with wrong   '{wrong_smi}': found {len(found)} results: {[r[0] for r in found]}")

            # Try BINARY comparison (force case-sensitive)
            sSql = f"SELECT compound_id FROM {db_name}.compound WHERE BINARY smiles_std = '{correct_smi}'"
            cur.execute(sSql)
            found_cs = cur.fetchall()
            print(f"  BINARY search  '{correct_smi}': found {len(found_cs)} results: {[r[0] for r in found_cs]}")

            sSql = f"SELECT compound_id FROM {db_name}.compound WHERE BINARY smiles_std = '{wrong_smi}'"
            cur.execute(sSql)
            found_cs = cur.fetchall()
            print(f"  BINARY search  '{wrong_smi}': found {len(found_cs)} results: {[r[0] for r in found_cs]}")
    except Exception as e:
        print(f"\n  {db_name}: {e}")


print()
print("=" * 100)
print("TEST 2: Registration timeline - when were these compounds first created?")
print("=" * 100)

for db_name in ['bcpvs', 'bcpvs_test']:
    try:
        for cid, info in sorted(AFFECTED_COMPOUNDS.items()):
            sSql = f"""SELECT compound_id, created_date, smiles_std
                       FROM {db_name}.compound
                       WHERE compound_id = '{cid}'"""
            cur.execute(sSql)
            rows = cur.fetchall()
            if rows:
                for (compound_id, created_date, smi) in rows:
                    print(f"  {db_name}: {compound_id} created={created_date} smiles_std={smi}")

            # Also check batches
            sSql = f"""SELECT compound_id, notebook_ref, submittal_date, chemreg_regno
                       FROM {db_name}.batch
                       WHERE compound_id = '{cid}'
                       ORDER BY submittal_date"""
            cur.execute(sSql)
            batches = cur.fetchall()
            for (compound_id, notebook_ref, submit_date, chemreg_regno) in batches:
                print(f"    Batch: {notebook_ref} submitted={submit_date} regno={chemreg_regno}")
    except Exception as e:
        print(f"  {db_name}: {e}")
    print()


print()
print("=" * 100)
print("TEST 3: Check for other Python/RDKit installations")
print("=" * 100)

# Check conda environments
try:
    result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True, timeout=10)
    print("Conda environments:")
    print(result.stdout)
except:
    print("  conda not available or not in PATH")

# Check which python the tornado server runs as
try:
    result = subprocess.run(['pgrep', '-af', 'angular.py'], capture_output=True, text=True, timeout=5)
    if result.stdout.strip():
        print("Running angular.py processes:")
        print(f"  {result.stdout.strip()}")
    else:
        print("No angular.py process currently running")
except:
    print("  Could not check running processes")

# Check which python the server uses (from any supervisor/systemd config)
for svc_file in ['/etc/supervisor/conf.d/scarab.conf',
                 '/etc/systemd/system/scarab.service',
                 os.path.join(parent_dir, 'start_server.sh'),
                 os.path.join(parent_dir, 'run.sh')]:
    if os.path.exists(svc_file):
        print(f"\n  Found service config: {svc_file}")
        with open(svc_file) as f:
            print(f"  Contents:\n  {f.read()}")


print()
print("=" * 100)
print("TEST 4: Check the JCMOL structure tables for the affected compounds")
print("        The molfile stored in JCMOL came from cleanStructureRDKit's Molcart output")
print("=" * 100)

for db_name in ['bcpvs', 'bcpvs_test']:
    try:
        for cid in ['CBK278468', 'CBK714175', 'CBK028998']:
            info = AFFECTED_COMPOUNDS[cid]
            # Read the molfile stored in JCMOL_MOLTABLE
            sSql = f"""SELECT compound_id, bin2smiles(mol, 'mol') as smi
                       FROM {db_name}.JCMOL_MOLTABLE_MOL
                       WHERE compound_id = '{cid}'"""
            cur.execute(sSql)
            rows = cur.fetchall()
            if rows:
                jcmol_smi = rows[0][1]
                if isinstance(jcmol_smi, bytes):
                    jcmol_smi = jcmol_smi.decode()
                print(f"  {db_name} JCMOL_MOLTABLE_MOL {cid}:")
                print(f"    SMILES from JCMOL structure: {jcmol_smi}")
                print(f"    compound.smiles_std:         {info['wrong_smi']}")
                print(f"    Correct SMILES:              {info['correct_smi']}")
                # Is the JCMOL structure also wrong?
                if jcmol_smi and info['correct_smi'] not in jcmol_smi:
                    print(f"    *** JCMOL structure is also WRONG ***")
            else:
                print(f"  {db_name} JCMOL_MOLTABLE_MOL: {cid} not found")
    except Exception as e:
        print(f"  {db_name}: {e}")


print()
print("=" * 100)
print("ANALYSIS")
print("=" * 100)
print("""
SCENARIO: The wrong smiles_std was stored by a PREVIOUS RDKit version.

Evidence:
1. Current RDKit 2022.3 produces CORRECT SMILES from these molfiles
2. The compound table has WRONG smiles_std
3. RDKit 2022.3 is deterministic (won't produce wrong results randomly)
4. The wrong forms all show aromatic <-> saturated ring swaps for amides
   - This is a well-known RDKit TautomerEnumerator bug from older versions

Timeline hypothesis:
1. File 1 registered with OLD buggy RDKit → produced wrong smiles_std
   → stored in compound table → matched existing wrong entries
2. RDKit upgraded to 2022.3 (or conda env changed)  
3. File 2 registered with NEW RDKit 2022.3 → produced correct smiles_std
   → correct SMILES doesn't match wrong stored value → registered as new

Key question: Was RDKit upgraded, or was the server conda environment
changed between the File 1 and File 2 registrations?

If the collation is case-insensitive, there's another twist:
   Even with wrong SMILES, the search might still match (c vs C same in CI).
   Check TEST 1 results carefully.
""")
