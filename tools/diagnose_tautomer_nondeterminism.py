#!/usr/bin/env python3
"""
Trace the cleanStructureRDKit pipeline step by step for the paired compounds
that got incorrectly merged. Also test for non-determinism across processes.
"""
import sys
import os
import subprocess
import tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

print(f"RDKit version: {Chem.rdBase.rdkitVersion}")

# Pairs of supplier compounds (different ext_ids) that got merged into the same compound_id
# Each pair: (compound_id, ext_id_1, smiles_1, ext_id_2, smiles_2)
MERGED_PAIRS = [
    ('CBK714175', 'Z32016367', 'CNC(=O)c1cccnc1', 'Z285662830', 'CNC(=O)C1CCCNC1.Cl'),
    ('CBK714200', 'Z31490287', 'O=C(NCC1CCCO1)c1cccc(Cl)c1', 'Z31407959', 'O=C(NCc1ccco1)c1cccc(Cl)c1'),
    ('CBK291635', 'Z27782625', 'O=C(Nc1ccccc1)c1cccnc1', 'Z31484247', 'O=C(NC1CCCCC1)c1cccnc1'),
    ('CBK699880', 'Z27749575_sdf', 'O=C(NCc1ccccc1)C1CCCCC1', None, None),  # The second is from a different source
]


def trace_pipeline(smiles_input, label):
    """Trace each step of the cleanStructureRDKit pipeline."""
    print(f"\n  [{label}] Input SMILES: {smiles_input}")
    
    mol = Chem.MolFromSmiles(smiles_input)
    if mol is None:
        print(f"  Parse failed!")
        return None
    
    raw = Chem.MolToSmiles(mol)
    print(f"  After parse:           {raw}")
    
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    
    clean_mol = rdMolStandardize.Cleanup(mol, params)
    print(f"  After Cleanup:         {Chem.MolToSmiles(clean_mol)}")
    
    parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    print(f"  After FragmentParent:  {Chem.MolToSmiles(parent_mol)}")
    
    uncharged = rdMolStandardize.Uncharger().uncharge(parent_mol)
    print(f"  After Uncharge:        {Chem.MolToSmiles(uncharged)}")
    
    te = rdMolStandardize.TautomerEnumerator(params)
    
    # Enumerate ALL tautomers
    tautomers = list(te.Enumerate(uncharged))
    print(f"  Tautomers enumerated:  {len(tautomers)}")
    for i, t in enumerate(tautomers):
        tsmi = Chem.MolToSmiles(t)
        print(f"    Tautomer {i}: {tsmi}")
    
    canon = te.Canonicalize(uncharged)
    final = Chem.MolToSmiles(canon)
    print(f"  After Canonicalize:    {final}")
    
    return final


print("=" * 100)
print("PART 1: Step-by-step pipeline trace for merged pairs")
print("=" * 100)

for cid, ext1, smi1, ext2, smi2 in MERGED_PAIRS:
    if ext2 is None:
        continue
    print(f"\n{'='*80}")
    print(f"{cid}: Two different supplier compounds merged")
    print(f"{'='*80}")
    
    result1 = trace_pipeline(smi1, f"{ext1}")
    result2 = trace_pipeline(smi2, f"{ext2}")
    
    if result1 and result2:
        print(f"\n  RESULT: {ext1}={result1}, {ext2}={result2}")
        print(f"  Same stdSMILES? {result1 == result2}")
        if result1 == result2:
            print(f"  *** BUG REPRODUCED: Two different molecules canonicalize to the same SMILES! ***")


print("\n\n")
print("=" * 100)
print("PART 2: Test non-determinism across process restarts (ASLR hypothesis)")
print("=" * 100)
print("Running cleanStructureRDKit in separate child processes...\n")

CHILD_SCRIPT = '''
import sys
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

smiles = sys.argv[1]
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("PARSE_FAILED")
    sys.exit(0)

params = rdMolStandardize.CleanupParameters()
params.tautomerRemoveSp3Stereo = False
params.tautomerRemoveBondStereo = False
params.tautomerRemoveIsotopicHs = False

clean_mol = rdMolStandardize.Cleanup(mol, params)
parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
uncharged = rdMolStandardize.Uncharger().uncharge(parent_mol)
te = rdMolStandardize.TautomerEnumerator(params)
canon = te.Canonicalize(uncharged)
print(Chem.MolToSmiles(canon))
'''

# Write child script to temp file
child_script_path = tempfile.mktemp(suffix='.py')
with open(child_script_path, 'w') as f:
    f.write(CHILD_SCRIPT)

test_smiles = [
    ('CNC(=O)C1CCCNC1.Cl', 'Piperidine amide + HCl'),
    ('CNC(=O)C1CCCNC1', 'Piperidine amide'),
    ('CNC(=O)c1ccncc1', 'Pyridine amide'),
    ('O=C(NCc1ccco1)c1cccc(Cl)c1', 'Furanylmethyl amide'),
    ('O=C(NCC1CCCO1)c1cccc(Cl)c1', 'THF-methyl amide'),
    ('O=C(NC1CCCCC1)c1cccnc1', 'Cyclohexylamine pyridine'),
    ('O=C(Nc1ccccc1)c1cccnc1', 'Aniline pyridine'),
    ('NC(=O)C1CCNCC1', 'Isonipecotamide'),
    ('NC(=O)c1ccncc1', 'Nicotinamide'),
    ('O=C(NC1CCCCC1)c1ccccc1', 'N-cyclohexylbenzamide'),
    ('O=C(Nc1ccccc1)c1ccccc1', 'Benzanilide'),
    ('O=C(NCC1CCCO1)Nc1ccccc1', 'THF-methyl urea'),
    ('O=C(NCC1CCCO1)NC1CCCCC1', 'THF-methyl cyclohexyl urea'),
]

NUM_TRIALS = 20

for smi, name in test_smiles:
    results = set()
    for trial in range(NUM_TRIALS):
        try:
            result = subprocess.check_output(
                [sys.executable, child_script_path, smi],
                stderr=subprocess.DEVNULL,
                timeout=30
            ).decode().strip()
            results.add(result)
        except Exception as e:
            results.add(f"ERROR: {e}")
    
    if len(results) == 1:
        only_result = list(results)[0]
        # Check if the canonical SMILES matches the input
        input_mol = Chem.MolFromSmiles(smi)
        params = rdMolStandardize.CleanupParameters()
        clean = rdMolStandardize.Cleanup(input_mol, params)
        parent = rdMolStandardize.FragmentParent(clean, params)
        uncharged_input = rdMolStandardize.Uncharger().uncharge(parent)
        input_canonical = Chem.MolToSmiles(uncharged_input)
        
        if only_result != input_canonical:
            print(f"  {name}: CONSISTENT but CHANGED: {input_canonical} → {only_result}")
        else:
            print(f"  {name}: deterministic → {only_result}")
    else:
        print(f"  {name}: *** NON-DETERMINISTIC! *** {len(results)} different results:")
        for r in sorted(results):
            print(f"    → {r}")

# Cleanup
os.unlink(child_script_path)


print("\n\n")
print("=" * 100)
print("PART 3: Test with MOLFILE input (as actual server does)")
print("=" * 100)

import MySQLdb
import config

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

# Get the actual molfiles from the DB for the affected regnos
AFFECTED_REGNOS = [
    (4052380, 'EH8331106', 'CBK714175', 'CNC(=O)C1CCCNC1.Cl'),  # Wrong (piperidine stored as pyridine)
    (4052279, 'EH8331005', 'CBK714175', 'CNC(=O)c1cccnc1'),      # Correct (pyridine stored correctly)
    (4052406, 'EH8331132', 'CBK714200', 'O=C(NCc1ccco1)c1cccc(Cl)c1'),  # Wrong
    (4052310, 'EH8331036', 'CBK714200', 'O=C(NCC1CCCO1)c1cccc(Cl)c1'),  # Correct
]

CHILD_MOLFILE_SCRIPT = '''
import sys
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

molfile = open(sys.argv[1]).read()
mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
if mol is None:
    print("PARSE_FAILED")
    sys.exit(0)

params = rdMolStandardize.CleanupParameters()
params.tautomerRemoveSp3Stereo = False
params.tautomerRemoveBondStereo = False
params.tautomerRemoveIsotopicHs = False

clean_mol = rdMolStandardize.Cleanup(mol, params)
parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
uncharged = rdMolStandardize.Uncharger().uncharge(parent_mol)
te = rdMolStandardize.TautomerEnumerator(params)
canon = te.Canonicalize(uncharged)
print(Chem.MolToSmiles(canon))
'''

child_script_path = tempfile.mktemp(suffix='.py')
with open(child_script_path, 'w') as f:
    f.write(CHILD_MOLFILE_SCRIPT)

for regno, jpage, cid, expected_smiles in AFFECTED_REGNOS:
    cur.execute("SELECT MOLFILE FROM chem_reg.chem_info WHERE regno = %s", (regno,))
    row = cur.fetchone()
    if not row or not row[0]:
        print(f"  {jpage}: No MOLFILE")
        continue
    
    molfile = row[0]
    
    # Write molfile to temp file
    molfile_path = tempfile.mktemp(suffix='.mol')
    with open(molfile_path, 'w') as f:
        f.write(molfile)
    
    results = set()
    for trial in range(NUM_TRIALS):
        try:
            result = subprocess.check_output(
                [sys.executable, child_script_path, molfile_path],
                stderr=subprocess.DEVNULL,
                timeout=30
            ).decode().strip()
            results.add(result)
        except Exception as e:
            results.add(f"ERROR: {e}")
    
    os.unlink(molfile_path)
    
    if len(results) == 1:
        only_result = list(results)[0]
        print(f"  {jpage} ({cid}): deterministic → {only_result}")
        if only_result != expected_smiles:
            # Check if it matches the wrong one
            print(f"    Expected raw: {expected_smiles}")
    else:
        print(f"  {jpage} ({cid}): *** NON-DETERMINISTIC! ***")
        for r in sorted(results):
            print(f"    → {r}")

os.unlink(child_script_path)
cur.close()
conn.close()
