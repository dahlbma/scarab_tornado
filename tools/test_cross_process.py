"""
Test to run ON THE SERVER to check if RDKit's TautomerEnumerator produces
different results across separate Python processes.

The C++ unordered containers (std::unordered_set, std::unordered_map)
used internally by TautomerEnumerator iterate in an order that depends on
memory addresses, which change with ASLR on each process start.

This script spawns multiple Python sub-processes, each running
cleanStructureRDKit on the same molfile, and checks if they produce
different SMILES.
"""
import subprocess
import sys
import os
import re
import io
import tempfile

AFFECTED_SUPPLIER_IDS = [
    'Z104477230', 'Z199507048', 'Z27749575', 'Z285662830',
    'Z31407959', 'Z31484244', 'Z31484247', 'Z31484290',
    'Z32016358', 'Z44585993',
]

# Known wrong values from the compound table
DB_REGISTERED_SMILES = {
    'Z104477230': 'NC(=O)c1ccncc1',
    'Z199507048': 'O=C(NC1CCNCC1)C1CCCCC1',
    'Z27749575':  'O=C(NCc1ccccc1)c1ccccc1',
    'Z285662830': 'CNC(=O)c1cccnc1',
    'Z31407959':  'O=C(NCC1CCCO1)c1cccc(Cl)c1',
    'Z31484244':  'O=C(Nc1ccccc1)c1cccc(O)c1',
    'Z31484247':  'O=C(Nc1ccccc1)c1cccnc1',
    'Z31484290':  'O=C(Nc1ccccc1)c1ccccc1',
    'Z32016358':  'CNC(=O)C1CCNCC1',
    'Z44585993':  'O=C(NCC1CCCO1)NC1CCCCC1',
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


# Child process script: reads a molfile from a temp file, runs cleanStructureRDKit,
# prints the SMILES
CHILD_SCRIPT = '''
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

try:
    clean_mol = rdMolStandardize.Cleanup(mol, params)
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(parent_clean_mol)
    te = rdMolStandardize.TautomerEnumerator(params)
    canonical = te.Canonicalize(uncharged)
    smiles = Chem.MolToSmiles(canonical)
    print(smiles)
except Exception as e:
    print(f"ERROR: {e}")
'''


def main():
    from rdkit import Chem
    print(f"RDKit version: {Chem.rdBase.rdkitVersion}")
    print(f"Python: {sys.executable}")
    print()

    # Find SDF files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    for candidate in [os.path.join(parent_dir, 'chemReg'), '.', '..', '../chemReg']:
        test_file = os.path.join(candidate, 'CBCS0523_CC04606_208cpds_Carlsson.sdf')
        if os.path.exists(test_file):
            sdf_dir = candidate
            break
    else:
        print("ERROR: Cannot find SDF files.")
        sys.exit(1)

    file1 = os.path.join(sdf_dir, 'CBCS0523_CC04606_208cpds_Carlsson.sdf')
    mols1 = parse_sdf(file1)

    # Write child script to temp file
    child_script_path = os.path.join(tempfile.gettempdir(), 'rdkit_tautomer_child.py')
    with open(child_script_path, 'w') as f:
        f.write(CHILD_SCRIPT)

    N_RUNS = 20
    print(f"Spawning {N_RUNS} separate Python processes for each molecule...")
    print(f"Each process has different ASLR (Address Space Layout Randomization)")
    print(f"which can affect C++ std::unordered_container iteration order.")
    print()

    print("=" * 100)
    print("CROSS-PROCESS NON-DETERMINISM TEST")
    print("=" * 100)

    any_non_deterministic = False

    for sid in AFFECTED_SUPPLIER_IDS:
        molblock = mols1.get(sid)
        if not molblock:
            continue

        # Write molfile to temp file
        mol_path = os.path.join(tempfile.gettempdir(), f'molfile_{sid}.mol')
        with open(mol_path, 'w') as f:
            f.write(molblock)

        results = {}
        for i in range(N_RUNS):
            try:
                result = subprocess.run(
                    [sys.executable, child_script_path, mol_path],
                    capture_output=True, text=True, timeout=30
                )
                smi = result.stdout.strip()
                if smi:
                    results[smi] = results.get(smi, 0) + 1
            except Exception as e:
                results[f"ERROR: {e}"] = results.get(f"ERROR: {e}", 0) + 1

        db_smi = DB_REGISTERED_SMILES.get(sid, '???')

        if len(results) > 1:
            any_non_deterministic = True
            print(f"\n  {sid}: *** NON-DETERMINISTIC across processes! ***")
            print(f"    DB registered SMILES: {db_smi}")
            for smi, count in sorted(results.items(), key=lambda x: -x[1]):
                marker = ""
                if smi == db_smi:
                    marker = " ← matches DB (the WRONG form)"
                print(f"    {count:3d}/{N_RUNS} runs: {smi}{marker}")
        else:
            smi = list(results.keys())[0] if results else None
            matches_db = " (matches DB)" if smi == db_smi else f" (DB has: {db_smi})"
            print(f"  {sid}: Deterministic across {N_RUNS} processes: {smi}{matches_db}")

        # Cleanup
        os.remove(mol_path)

    os.remove(child_script_path)

    print()
    print("=" * 100)
    if any_non_deterministic:
        print("CONCLUSION: TautomerEnumerator IS non-deterministic across process restarts!")
        print("This is caused by ASLR changing memory addresses, which affects")
        print("C++ std::unordered_container iteration order in RDKit's tautomer scoring.")
        print()
        print("At registration time, the server process happened to produce the wrong")
        print("aromatic/saturated SMILES. After server restart, it produces correct SMILES.")
        print("The second file registration (after restart) couldn't find the wrong SMILES")
        print("in the compound table → molecules registered as new duplicates.")
    else:
        print("CONCLUSION: TautomerEnumerator appears deterministic across processes")
        print("with current RDKit version. The bug may have been triggered by a")
        print("different condition at registration time (different RDKit version,")
        print("Molcart version, or other environmental factor).")
    print("=" * 100)


if __name__ == '__main__':
    main()
