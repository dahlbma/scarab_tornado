"""
Check when the tornado server process was started and what RDKit
version exists in the other anaconda installation.

Something changed between 2026-03-18 (File 1, wrong SMILES produced)
and now (correct SMILES). The server says RDKit wasn't updated, but:
- /home/scarab_data/anaconda3 exists as another Python installation
- The server process might have been restarted with a different env
"""
import subprocess
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, parent_dir)

from rdkit import Chem
print(f"Current RDKit: {Chem.rdBase.rdkitVersion}")
print(f"Current Python: {sys.executable}")
print()

# 1. When was the tornado server process started?
print("=" * 100)
print("TEST 1: When was the tornado server process started?")
print("=" * 100)
try:
    result = subprocess.run(['pgrep', '-af', 'angular.py'], capture_output=True, text=True, timeout=5)
    if result.stdout.strip():
        for line in result.stdout.strip().split('\n'):
            pid = line.split()[0]
            print(f"  PID: {pid}")
            print(f"  Command: {line}")
            # Get process start time
            result2 = subprocess.run(['ps', '-o', 'lstart=', '-p', pid],
                                     capture_output=True, text=True, timeout=5)
            print(f"  Started: {result2.stdout.strip()}")
            # Get process elapsed time
            result3 = subprocess.run(['ps', '-o', 'etime=', '-p', pid],
                                     capture_output=True, text=True, timeout=5)
            print(f"  Elapsed: {result3.stdout.strip()}")
    else:
        print("  No angular.py process running")
except Exception as e:
    print(f"  Error: {e}")


# 2. Check RDKit in /home/scarab_data/anaconda3
print()
print("=" * 100)
print("TEST 2: Check RDKit in /home/scarab_data/anaconda3")
print("=" * 100)

other_python_candidates = [
    '/home/scarab_data/anaconda3/bin/python',
    '/home/scarab_data/anaconda3/bin/python3',
]

for py_path in other_python_candidates:
    if os.path.exists(py_path):
        print(f"\n  Found: {py_path}")
        try:
            result = subprocess.run(
                [py_path, '-c',
                 'from rdkit import Chem; print(f"RDKit: {Chem.rdBase.rdkitVersion}")'],
                capture_output=True, text=True, timeout=10
            )
            if result.stdout.strip():
                print(f"  {result.stdout.strip()}")
            if result.stderr.strip():
                print(f"  stderr: {result.stderr.strip()[:200]}")
        except Exception as e:
            print(f"  Error: {e}")
    else:
        print(f"  Not found: {py_path}")

# Also check conda envs in that installation
other_conda = '/home/scarab_data/anaconda3/bin/conda'
if os.path.exists(other_conda):
    print(f"\n  Conda envs in /home/scarab_data/anaconda3:")
    try:
        result = subprocess.run([other_conda, 'env', 'list'],
                                capture_output=True, text=True, timeout=10)
        print(f"  {result.stdout}")
    except Exception as e:
        print(f"  Error: {e}")

    # Check rdkit package specifically
    print(f"  RDKit package in /home/scarab_data/anaconda3:")
    try:
        result = subprocess.run([other_conda, 'list', 'rdkit'],
                                capture_output=True, text=True, timeout=10)
        print(f"  {result.stdout}")
    except Exception as e:
        print(f"  Error: {e}")


# 3. Check conda revision history for my-rdkit-env
print()
print("=" * 100)
print("TEST 3: Conda revision history for my-rdkit-env")
print("        Shows all package installs/updates with dates")
print("=" * 100)
try:
    result = subprocess.run(['conda', 'list', '--revisions', '-n', 'my-rdkit-env'],
                            capture_output=True, text=True, timeout=30)
    # Only show the last few revisions and any rdkit-related changes
    lines = result.stdout.strip().split('\n')
    print(f"  Total lines: {len(lines)}")
    # Print last 60 lines
    for line in lines[-60:]:
        print(f"  {line}")
except Exception as e:
    print(f"  Error: {e}")


# 4. Check when rdkit .so files were last modified
print()
print("=" * 100)
print("TEST 4: RDKit shared library modification dates")
print("=" * 100)
try:
    rdkit_path = os.path.dirname(os.path.dirname(Chem.__file__))
    print(f"  RDKit installed at: {rdkit_path}")

    # Find .so files
    result = subprocess.run(
        ['find', rdkit_path, '-name', '*.so', '-newer', '/etc/hostname', '-ls'],
        capture_output=True, text=True, timeout=10
    )
    # Just show timestamps of a few key files
    result2 = subprocess.run(
        ['ls', '-la', os.path.join(os.path.dirname(Chem.__file__), 'rdBase.so')],
        capture_output=True, text=True, timeout=5
    )
    if result2.stdout.strip():
        print(f"  {result2.stdout.strip()}")

    # Check Chem module
    result3 = subprocess.run(
        ['ls', '-la', os.path.join(os.path.dirname(Chem.__file__), 'Chem', 'rdchem.so')],
        capture_output=True, text=True, timeout=5
    )
    if result3.stdout.strip():
        print(f"  {result3.stdout.strip()}")

    # Check MolStandardize
    msdir = os.path.join(os.path.dirname(Chem.__file__), 'Chem', 'MolStandardize')
    result4 = subprocess.run(
        ['ls', '-la', os.path.join(msdir, 'rdMolStandardize.so')],
        capture_output=True, text=True, timeout=5
    )
    if result4.stdout.strip():
        print(f"  {result4.stdout.strip()}")
except Exception as e:
    print(f"  Error: {e}")


# 5. Check if the running server process uses the same RDKit
print()
print("=" * 100)
print("TEST 5: What RDKit .so is loaded by the RUNNING tornado server?")
print("=" * 100)
try:
    result = subprocess.run(['pgrep', '-f', 'angular.py'], capture_output=True, text=True, timeout=5)
    if result.stdout.strip():
        pid = result.stdout.strip().split('\n')[0]
        # Check /proc/PID/maps for rdkit shared libraries
        maps_file = f'/proc/{pid}/maps'
        if os.path.exists(maps_file):
            result2 = subprocess.run(['grep', '-i', 'rdkit', maps_file],
                                     capture_output=True, text=True, timeout=5)
            if result2.stdout.strip():
                # Deduplicate the library paths
                libs = set()
                for line in result2.stdout.strip().split('\n'):
                    parts = line.split()
                    if len(parts) >= 6:
                        libs.add(parts[-1])
                print(f"  RDKit libraries loaded by server (PID {pid}):")
                for lib in sorted(libs)[:10]:
                    print(f"    {lib}")
                    # Show file date
                    result3 = subprocess.run(['ls', '-la', lib],
                                             capture_output=True, text=True, timeout=5)
                    if result3.stdout.strip():
                        print(f"      {result3.stdout.strip()}")
            else:
                print(f"  No RDKit libraries found in /proc/{pid}/maps")
        else:
            print(f"  Cannot access {maps_file}")
    else:
        print("  No angular.py process running")
except Exception as e:
    print(f"  Error: {e}")


# 6. When was the Molcart/scarab container last started?
print()
print("=" * 100)
print("TEST 6: Molcart/Scarab Podman container status")
print("=" * 100)
try:
    result = subprocess.run(['podman', 'ps', '-a', '--format',
                             '{{.Names}} {{.Image}} {{.Created}} {{.Status}}'],
                            capture_output=True, text=True, timeout=10)
    if result.stdout.strip():
        print(f"  {result.stdout.strip()}")
    else:
        print("  No podman containers found (might need root)")
        # Try with sudo or check systemctl
        result2 = subprocess.run(['systemctl', 'status', 'scarab.service'],
                                 capture_output=True, text=True, timeout=10)
        if result2.stdout.strip():
            for line in result2.stdout.strip().split('\n')[:15]:
                print(f"  {line}")
except Exception as e:
    print(f"  Error: {e}")
