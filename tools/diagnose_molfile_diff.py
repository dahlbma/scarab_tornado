#!/usr/bin/env python3
"""
Compare the MOLFILE stored in chem_reg.chem_info with the original SDF molblock.
This will show if the MOLFILE was corrupted during storage.
"""
import sys
import os
import re
import io
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import MySQLdb
import config
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

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


def parse_sdf(filename):
    """Parse SDF file - return ALL entries (not deduped by external_id if multiple exist)."""
    results = []
    with open(filename, 'rb') as f:
        content = f.read()
    f = io.BytesIO(content)
    mol_index = 0
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
        # Extract external ID
        pattern = r'>\s*<External ID>\s*\n(\S+)'
        match = re.search(pattern, mol_str)
        ext_id = match.group(1).strip() if match else f"UNKNOWN_{mol_index}"
        # Extract molblock
        m_end_idx = mol_str.find('M  END')
        molblock = mol_str[:m_end_idx + 6] if m_end_idx >= 0 else None
        results.append({'index': mol_index, 'ext_id': ext_id, 'molblock': molblock, 'full_record': mol_str})
        mol_index += 1
    return results


# The affected regnos and their batch info
AFFECTED_REGNOS = [
    (4052279, 'EH8331005', 'CBK714175'),  # Correct
    (4052380, 'EH8331106', 'CBK714175'),  # Wrong
    (4052310, 'EH8331036', 'CBK714200'),  # Correct  
    (4052406, 'EH8331132', 'CBK714200'),  # Wrong
    (4052333, 'EH8331059', 'CBK291635'),  # Wrong
    (4052346, 'EH8331072', 'CBK291635'),  # Correct
    (4052300, 'EH8331026', 'CBK699880'),  # Wrong
    (4052394, 'EH8331120', 'CBK028998'),  # Wrong
    (4052400, 'EH8331126', 'CBK228517'),  # Wrong
    (4052285, 'EH8331011', 'CBK278635'),  # Wrong
]

print("=" * 120)
print("PART 1: Get external_id for each affected batch from chem_reg")
print("=" * 120)

for regno, jpage, cid in AFFECTED_REGNOS:
    cur.execute("SELECT EXTERNAL_ID, COMPOUND_ID, sdfile_sequence FROM chem_reg.chem_info WHERE regno = %s", (regno,))
    row = cur.fetchone()
    if row:
        print(f"  {jpage} (regno={regno}) → compound={row[1]}, ext_id={row[0]}, sdf_seq={row[2]}")


print("\n\n")
print("=" * 120)
print("PART 2: Compare stored MOLFILE in DB vs original SDF for affected regnos")
print("=" * 120)

sdf_dir = '../chemReg'
file1 = f'{sdf_dir}/CBCS0523_CC04606_208cpds_Carlsson.sdf'
file2 = f'{sdf_dir}/CBCS0523_CC04640_192cpds_Carlsson.sdf'

print("Parsing SDF files...")
mols1 = parse_sdf(file1)
mols2 = parse_sdf(file2)

# Build lookup by external_id (keeping all occurrences)
sdf_by_extid = {}
for m in mols1:
    sdf_by_extid.setdefault(m['ext_id'], []).append(('File1', m))
for m in mols2:
    sdf_by_extid.setdefault(m['ext_id'], []).append(('File2', m))

print(f"File 1: {len(mols1)} molecules")
print(f"File 2: {len(mols2)} molecules")
print()

for regno, jpage, cid in AFFECTED_REGNOS:
    print(f"\n--- {jpage} (regno={regno}, compound={cid}) ---")
    
    # Get stored molfile from DB
    cur.execute("SELECT MOLFILE, EXTERNAL_ID FROM chem_reg.chem_info WHERE regno = %s", (regno,))
    row = cur.fetchone()
    if not row or not row[0]:
        print("  No MOLFILE in chem_reg!")
        continue
    db_molfile = row[0]
    ext_id = row[1]
    print(f"  External ID: {ext_id}")
    
    # Extract just the molblock from DB molfile (up to M  END)
    m_end_idx = db_molfile.find('M  END')
    if m_end_idx >= 0:
        db_molblock = db_molfile[:m_end_idx + 6]
    else:
        db_molblock = db_molfile
        print("  WARNING: No 'M  END' found in stored MOLFILE!")
    
    # RDKit SMILES from DB molfile
    mol_db = Chem.MolFromMolBlock(db_molblock, removeHs=False, sanitize=True)
    db_smiles = Chem.MolToSmiles(mol_db) if mol_db else "PARSE_FAILED"
    print(f"  DB MOLFILE → raw SMILES: {db_smiles}")
    
    # Find in SDF files
    if ext_id in sdf_by_extid:
        for source, sdf_entry in sdf_by_extid[ext_id]:
            sdf_molblock = sdf_entry['molblock']
            mol_sdf = Chem.MolFromMolBlock(sdf_molblock, removeHs=False, sanitize=True)
            sdf_smiles = Chem.MolToSmiles(mol_sdf) if mol_sdf else "PARSE_FAILED"
            
            # Compare molblocks
            db_lines = db_molblock.strip().split('\n')
            sdf_lines = sdf_molblock.strip().split('\n')
            
            identical = True
            for i, (dl, sl) in enumerate(zip(db_lines, sdf_lines)):
                if dl.rstrip() != sl.rstrip():
                    identical = False
                    if i > 0:  # Skip header line
                        print(f"  DIFF at line {i}:")
                        print(f"    DB:  '{dl.rstrip()}'")
                        print(f"    SDF: '{sl.rstrip()}'")
            
            if len(db_lines) != len(sdf_lines):
                identical = False
                print(f"  Different line count: DB={len(db_lines)}, SDF={len(sdf_lines)}")
            
            match_str = "IDENTICAL" if identical else "DIFFERENT!"
            smiles_match = "same" if db_smiles == sdf_smiles else "DIFFERENT!"
            print(f"  {source} (index={sdf_entry['index']}): molblock={match_str}, SMILES: DB={db_smiles} SDF={sdf_smiles} ({smiles_match})")
    else:
        print(f"  External ID {ext_id} not found in either SDF file")


print("\n\n")
print("=" * 120)
print("PART 3: Check for duplicate external IDs within each SDF file")
print("=" * 120)

from collections import Counter
ext_ids1 = [m['ext_id'] for m in mols1]
ext_ids2 = [m['ext_id'] for m in mols2]
dups1 = {k: v for k, v in Counter(ext_ids1).items() if v > 1}
dups2 = {k: v for k, v in Counter(ext_ids2).items() if v > 1}
common_ids = set(ext_ids1) & set(ext_ids2)
print(f"  Duplicate ext_ids in File 1: {dups1 if dups1 else 'None'}")
print(f"  Duplicate ext_ids in File 2: {dups2 if dups2 else 'None'}")
print(f"  Common ext_ids between files: {len(common_ids)}")
print(f"  File 1 unique: {len(set(ext_ids1) - common_ids)}")
print(f"  File 2 unique: {len(set(ext_ids2) - common_ids)}")


print("\n\n")
print("=" * 120)
print("PART 4: Look at DB molfile character analysis for a wrong batch")
print("=" * 120)

# Check for characters that could cause issues: backslashes, non-printable, etc.
for regno, jpage, cid in [(4052380, 'EH8331106', 'CBK714175'), (4052279, 'EH8331005', 'CBK714175')]:
    print(f"\n--- {jpage} (regno={regno}) ---")
    cur.execute("SELECT MOLFILE FROM chem_reg.chem_info WHERE regno = %s", (regno,))
    row = cur.fetchone()
    if not row or not row[0]:
        continue
    db_molfile = row[0]
    
    # Character analysis
    char_counts = {}
    for ch in db_molfile:
        if ord(ch) < 32 and ch not in '\n\r\t':
            key = f"ctrl-{ord(ch)}"
            char_counts[key] = char_counts.get(key, 0) + 1
        elif ord(ch) > 127:
            key = f"U+{ord(ch):04X}"
            char_counts[key] = char_counts.get(key, 0) + 1
        elif ch == '\\':
            char_counts['backslash'] = char_counts.get('backslash', 0) + 1
    if char_counts:
        print(f"  Special chars: {char_counts}")
    else:
        print(f"  No special characters found (only ASCII printable + newline/tab)")
    
    print(f"  Total length: {len(db_molfile)}")
    print(f"  First 5 lines:")
    for i, line in enumerate(db_molfile.split('\n')[:5]):
        print(f"    [{i}] '{line}'")

cur.close()
conn.close()
