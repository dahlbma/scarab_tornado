"""
Test script to compare molblocks between the two SDF files for the 10 affected 
supplier IDs, and trace through the RDKit standardization pipeline to identify
why duplicate detection failed.
"""
import re
import sys
sys.path.insert(0, '..')

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# The 10 affected supplier IDs
AFFECTED_SUPPLIER_IDS = [
    'Z104477230',
    'Z199507048',
    'Z27749575',
    'Z285662830',
    'Z31407959',
    'Z31484244',
    'Z31484247',
    'Z31484290',
    'Z32016358',
    'Z44585993',
]

def parse_sdf(filename):
    """Parse an SDF file and return a dict of {external_id: (molblock, tags)}.
    
    Simulates sdfreg.py's getNextMolecule exactly: reads line by line,
    replaces \\r\\n with \\n, removes single quotes, prepends 'id' to first line.
    """
    import io
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
        
        # Extract external ID from tags
        pattern = r'>\s*<External ID>\s*\n(\S+)'
        match = re.search(pattern, mol_str)
        if match:
            ext_id = match.group(1).strip()
        else:
            continue
        
        # Extract molblock (everything up to M  END)
        m_end_idx = mol_str.find('M  END')
        if m_end_idx >= 0:
            molblock = mol_str[:m_end_idx + 6]
        else:
            continue
        
        results[ext_id] = (molblock, mol_str)
    
    return results


def cleanStructureRDKit_local(molfile):
    """
    Local version of cleanStructureRDKit from the server's dbInterface.py.
    Reproduces the same standardization pipeline without the SQL calls.
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
    except Exception as e:
        print(f"  Failed Cleanup: {e}")
        return None, None
    
    try:
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    except Exception as e:
        print(f"  Failed FragmentParent: {e}")
        return None, None
    
    uncharger = rdMolStandardize.Uncharger()
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
    
    te = rdMolStandardize.TautomerEnumerator(params)
    taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)
    
    stdSMILES = Chem.MolToSmiles(taut_uncharged_parent_clean_mol)
    molblock_out = Chem.MolToMolBlock(taut_uncharged_parent_clean_mol)
    
    return molblock_out, stdSMILES


def smiles_without_tautomer(molfile):
    """
    Same pipeline as cleanStructureRDKit but WITHOUT the TautomerEnumerator step.
    """
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return None
    
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    
    try:
        clean_mol = rdMolStandardize.Cleanup(mol, params)
    except:
        return None
    
    try:
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    except:
        return None
    
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(parent_clean_mol)
    
    return Chem.MolToSmiles(uncharged)


def enumerate_tautomers(molfile):
    """Enumerate all tautomers for a molecule to see what the enumerator considers."""
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return []
    
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    
    clean_mol = rdMolStandardize.Cleanup(mol, params)
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(parent_clean_mol)
    
    te = rdMolStandardize.TautomerEnumerator(params)
    tautomers = te.Enumerate(uncharged)
    
    return [Chem.MolToSmiles(t) for t in tautomers]


def main():
    sdf_dir = '../chemReg'
    file1 = f'{sdf_dir}/CBCS0523_CC04606_208cpds_Carlsson.sdf'
    file2 = f'{sdf_dir}/CBCS0523_CC04640_192cpds_Carlsson.sdf'
    
    print("Parsing SDF files...")
    mols1 = parse_sdf(file1)
    mols2 = parse_sdf(file2)
    
    print(f"File 1: {len(mols1)} molecules")
    print(f"File 2: {len(mols2)} molecules")
    print()
    
    print("=" * 80)
    print("STEP 1: Compare molblocks between files for affected supplier IDs")
    print("=" * 80)
    
    for sid in AFFECTED_SUPPLIER_IDS:
        in_file1 = sid in mols1
        in_file2 = sid in mols2
        
        if in_file1 and in_file2:
            molblock1 = mols1[sid][0]
            molblock2 = mols2[sid][0]
            
            # Compare ignoring the first line (header with date)
            lines1 = molblock1.strip().split('\n')
            lines2 = molblock2.strip().split('\n')
            
            # Skip first line (molecule name) and compare rest
            atoms_bonds1 = '\n'.join(lines1[1:])
            atoms_bonds2 = '\n'.join(lines2[1:])
            
            if atoms_bonds1 == atoms_bonds2:
                status = "IDENTICAL (excluding header)"
            else:
                status = "DIFFERENT!"
                # Find which lines differ
                for i, (l1, l2) in enumerate(zip(lines1, lines2)):
                    if l1 != l2:
                        print(f"  Line {i}: File1='{l1}' File2='{l2}'")
            
            print(f"  {sid}: {status}")
        elif in_file1:
            print(f"  {sid}: Only in file 1")
        elif in_file2:
            print(f"  {sid}: Only in file 2")
        else:
            print(f"  {sid}: NOT FOUND in either file!")
    
    print()
    print("=" * 80)
    print("STEP 2: RDKit standardization pipeline comparison")
    print("=" * 80)
    
    for sid in AFFECTED_SUPPLIER_IDS:
        if sid not in mols1 or sid not in mols2:
            continue
        
        molblock1 = mols1[sid][0]
        molblock2 = mols2[sid][0]
        
        print(f"\n--- {sid} ---")
        
        # Direct SMILES from molblock (no standardization)
        mol1 = Chem.MolFromMolBlock(molblock1, removeHs=False, sanitize=True)
        mol2 = Chem.MolFromMolBlock(molblock2, removeHs=False, sanitize=True)
        raw_smi1 = Chem.MolToSmiles(mol1) if mol1 else "PARSE_FAILED"
        raw_smi2 = Chem.MolToSmiles(mol2) if mol2 else "PARSE_FAILED"
        
        print(f"  Raw SMILES (no standardization):")
        print(f"    File1: {raw_smi1}")
        print(f"    File2: {raw_smi2}")
        print(f"    Match: {raw_smi1 == raw_smi2}")
        
        # SMILES without tautomer step
        no_taut1 = smiles_without_tautomer(molblock1)
        no_taut2 = smiles_without_tautomer(molblock2)
        print(f"  SMILES (cleanup+uncharge, NO tautomer):")
        print(f"    File1: {no_taut1}")
        print(f"    File2: {no_taut2}")
        print(f"    Match: {no_taut1 == no_taut2}")
        
        # Full standardization (same as server)
        _, std_smi1 = cleanStructureRDKit_local(molblock1)
        _, std_smi2 = cleanStructureRDKit_local(molblock2)
        print(f"  SMILES (full standardization WITH tautomer):")
        print(f"    File1: {std_smi1}")
        print(f"    File2: {std_smi2}")
        print(f"    Match: {std_smi1 == std_smi2}")
        
        if raw_smi1 != std_smi1:
            print(f"  *** TAUTOMER CHANGED STRUCTURE: {raw_smi1} -> {std_smi1}")

    print()
    print("=" * 80)
    print("STEP 3: Tautomer enumeration for affected molecules")
    print("=" * 80)
    
    for sid in AFFECTED_SUPPLIER_IDS:
        if sid not in mols1:
            continue
        
        molblock1 = mols1[sid][0]
        print(f"\n--- {sid} ---")
        
        tautomers = enumerate_tautomers(molblock1)
        print(f"  Number of tautomers: {len(tautomers)}")
        for i, t in enumerate(tautomers):
            print(f"    Tautomer {i}: {t}")

    print()
    print("=" * 80)
    print("STEP 4: Check if getNextMolecule preprocessing affects the result")
    print("=" * 80)
    print("Testing with 'id' prepended to first line (as sdfreg.py does)...")
    
    for sid in AFFECTED_SUPPLIER_IDS[:3]:  # Test just a few
        if sid not in mols1:
            continue
        
        molblock1 = mols1[sid][0]
        
        # Simulate getNextMolecule's modifications
        lines = molblock1.split('\n')
        lines[0] = 'id' + lines[0]
        modified_molblock = '\n'.join(lines)
        
        # Also remove single quotes (as sdfreg does)
        modified_molblock = modified_molblock.replace("'", "")
        
        _, std_original = cleanStructureRDKit_local(molblock1)
        _, std_modified = cleanStructureRDKit_local(modified_molblock)
        
        print(f"  {sid}:")
        print(f"    Original: {std_original}")
        print(f"    Modified (id prefix, no quotes): {std_modified}")
        print(f"    Match: {std_original == std_modified}")


if __name__ == '__main__':
    main()
