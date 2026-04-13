#!/usr/bin/env python3
"""
Test: Determine the PROVENANCE of wrong smiles_std values.

For each affected compound, we check:
1. compound.smiles_std (what's stored)
2. RDKit SMILES from chem_info.molfile (original supplier molfile)
3. RDKit SMILES from JCMOL_MOLTABLE.mol (Molcart-generated molfile)
4. Molcart bin2smiles from JCMOL_MOLTABLE.mol
5. Whether smiles_std matches any of the above sources

This will reveal HOW the wrong smiles_std was set:
- If smiles_std matches RDKit-from-JCMOL → set by re-processing Molcart's wrong molfiles
- If smiles_std matches Molcart-bin2smiles → set by Molcart directly (migration script)
- If smiles_std matches RDKit-from-chem_info → set correctly during registration (shouldn't be wrong)
"""

import mysql.connector
import sys
sys.path.insert(0, '..')

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def rdkit_standardize(molfile):
    """Reproduce cleanStructureRDKit's SMILES generation."""
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return None
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
        canon = te.Canonicalize(uncharged)
        return Chem.MolToSmiles(canon)
    except Exception as e:
        return f"ERROR: {e}"


def main():
    db = mysql.connector.connect(
        host="127.0.0.1",
        user="mats",
        passwd="",
        database="bcpvs"
    )
    cur = db.cursor()

    affected_compounds = [
        'CBK028998', 'CBK228517', 'CBK278468', 'CBK278635',
        'CBK291635', 'CBK291792', 'CBK389719', 'CBK699880',
        'CBK714175', 'CBK714200'
    ]

    # Also get the "duplicate" compound IDs from File 2
    duplicate_compounds = [
        'CBK714368', 'CBK714369', 'CBK714372', 'CBK714370',
        'CBK714371', 'CBK714373', 'CBK714374', 'CBK714375',
        'CBK714376', 'CBK714377'
    ]

    print("=" * 120)
    print("TEST 1: Provenance of smiles_std for affected compounds")
    print("=" * 120)

    for cid in affected_compounds:
        print(f"\n{'='*80}")
        print(f"COMPOUND: {cid}")
        print(f"{'='*80}")

        # 1. Get smiles_std from compound table
        cur.execute(f"SELECT smiles_std, created_date FROM bcpvs.compound WHERE compound_id = '{cid}'")
        rows = cur.fetchall()
        if not rows:
            print(f"  NOT FOUND in compound table!")
            continue
        stored_smiles = rows[0][0]
        created_date = rows[0][1]
        print(f"  compound.smiles_std = {stored_smiles}")
        print(f"  created_date        = {created_date}")

        # 2. Get JCMOL molfile
        cur.execute(f"SELECT mol FROM bcpvs.JCMOL_MOLTABLE WHERE compound_id = '{cid}'")
        jcmol_rows = cur.fetchall()
        jcmol_molfile = jcmol_rows[0][0] if jcmol_rows else None

        # 3. Get earliest chem_info molfile for this compound
        cur.execute(f"""SELECT ci.molfile, ci.regno, ci.rdate, ci.jpage
                       FROM chem_reg.chem_info ci
                       WHERE ci.compound_id = '{cid}' AND ci.molfile IS NOT NULL
                       ORDER BY ci.regno ASC LIMIT 1""")
        ci_rows = cur.fetchall()
        ci_molfile = ci_rows[0][0] if ci_rows else None
        ci_regno = ci_rows[0][1] if ci_rows else None
        ci_rdate = ci_rows[0][2] if ci_rows else None
        ci_jpage = ci_rows[0][3] if ci_rows else None
        print(f"  Earliest chem_info: regno={ci_regno}, rdate={ci_rdate}, jpage={ci_jpage}")

        # 4. RDKit SMILES from chem_info molfile (original supplier)
        rdkit_from_ci = rdkit_standardize(ci_molfile) if ci_molfile else None
        print(f"  RDKit from chem_info = {rdkit_from_ci}")

        # 5. RDKit SMILES from JCMOL molfile (Molcart-generated)
        rdkit_from_jcmol = rdkit_standardize(jcmol_molfile) if jcmol_molfile else None
        print(f"  RDKit from JCMOL     = {rdkit_from_jcmol}")

        # 6. Molcart bin2smiles from JCMOL molfile
        if jcmol_molfile:
            try:
                esc_mol = jcmol_molfile.replace("'", "\\'").replace("\\", "\\\\")
                cur.execute(f"SELECT bin2smiles(mol2bin('{esc_mol}', 'mol'), 'mol')")
                mc_rows = cur.fetchall()
                molcart_smiles = mc_rows[0][0].decode() if mc_rows and mc_rows[0][0] else None
            except Exception as e:
                molcart_smiles = f"ERROR: {e}"
        else:
            molcart_smiles = None
        print(f"  Molcart bin2smiles   = {molcart_smiles}")

        # 7. Molcart bin2smiles using mol2bin from the stored smiles_std
        if stored_smiles:
            try:
                cur.execute(f"SELECT bin2smiles(mol2bin('{stored_smiles}', 'smiles'), 'mol')")
                mc2_rows = cur.fetchall()
                molcart_roundtrip = mc2_rows[0][0].decode() if mc2_rows and mc2_rows[0][0] else None
            except Exception as e:
                molcart_roundtrip = f"ERROR: {e}"
        else:
            molcart_roundtrip = None
        print(f"  Molcart roundtrip    = {molcart_roundtrip}")

        # 8. Compare
        print(f"\n  --- MATCHES ---")
        if stored_smiles and rdkit_from_ci:
            print(f"  stored == RDKit(chem_info)?  {stored_smiles == rdkit_from_ci}")
        if stored_smiles and rdkit_from_jcmol:
            print(f"  stored == RDKit(JCMOL)?      {stored_smiles == rdkit_from_jcmol}")
        if stored_smiles and molcart_smiles and not molcart_smiles.startswith("ERROR"):
            print(f"  stored == Molcart(JCMOL)?    {stored_smiles == molcart_smiles}")

    print("\n\n")
    print("=" * 120)
    print("TEST 2: Check if ALL old compounds (created_date IS NULL) have smiles_std")
    print("=" * 120)

    cur.execute("""SELECT
        COUNT(*) as total,
        SUM(CASE WHEN smiles_std IS NOT NULL AND smiles_std != '' THEN 1 ELSE 0 END) as has_smiles,
        SUM(CASE WHEN smiles_std IS NULL OR smiles_std = '' THEN 1 ELSE 0 END) as no_smiles
        FROM bcpvs.compound
        WHERE created_date IS NULL""")
    row = cur.fetchall()[0]
    print(f"  Old compounds (created_date IS NULL): total={row[0]}, has_smiles={row[1]}, no_smiles={row[2]}")

    cur.execute("""SELECT
        COUNT(*) as total,
        SUM(CASE WHEN smiles_std IS NOT NULL AND smiles_std != '' THEN 1 ELSE 0 END) as has_smiles,
        SUM(CASE WHEN smiles_std IS NULL OR smiles_std = '' THEN 1 ELSE 0 END) as no_smiles
        FROM bcpvs.compound
        WHERE created_date IS NOT NULL""")
    row = cur.fetchall()[0]
    print(f"  New compounds (created_date IS NOT NULL): total={row[0]}, has_smiles={row[1]}, no_smiles={row[2]}")

    print("\n\n")
    print("=" * 120)
    print("TEST 3: Sample old compounds - do their smiles_std match RDKit from JCMOL?")
    print("=" * 120)

    cur.execute("""SELECT c.compound_id, c.smiles_std, j.mol
                   FROM bcpvs.compound c
                   JOIN bcpvs.JCMOL_MOLTABLE j ON c.compound_id = j.compound_id
                   WHERE c.created_date IS NULL
                   AND c.smiles_std IS NOT NULL AND c.smiles_std != ''
                   LIMIT 20""")
    sample_rows = cur.fetchall()

    match_jcmol = 0
    mismatch_jcmol = 0
    for row in sample_rows:
        cid, stored, jmol = row
        rdkit_jcmol = rdkit_standardize(jmol) if jmol else None
        matches = (stored == rdkit_jcmol)
        if matches:
            match_jcmol += 1
        else:
            mismatch_jcmol += 1
            print(f"  {cid}: stored={stored}, RDKit(JCMOL)={rdkit_jcmol}")
    print(f"  Results: {match_jcmol} match, {mismatch_jcmol} mismatch out of {len(sample_rows)} sampled")


    print("\n\n")
    print("=" * 120)
    print("TEST 4: For File 2 duplicate compounds - what smiles_std was assigned?")
    print("=" * 120)
    
    for cid in duplicate_compounds:
        cur.execute(f"SELECT smiles_std, created_date FROM bcpvs.compound WHERE compound_id = '{cid}'")
        rows = cur.fetchall()
        if rows:
            print(f"  {cid}: smiles_std={rows[0][0]}, created={rows[0][1]}")

    print("\n\n")
    print("=" * 120)
    print("TEST 5: Molcart mol2bin/bin2mol SMILES transformation test")
    print("  Does Molcart change SMILES when doing mol2bin(smiles,'smiles') -> bin2mol?")
    print("=" * 120)

    test_smiles = [
        'NC(=O)C1CCNCC1',       # Correct piperidine
        'NC(=O)c1ccncc1',       # Wrong pyridine
        'CNC(=O)C1CCCNC1',     # Correct piperidine
        'CNC(=O)c1cccnc1',     # Wrong pyridine
        'O=C(NCc1ccccc1)C1CCCCC1',  # Correct cyclohexane
        'O=C(NCc1ccccc1)c1ccccc1',  # Wrong benzene
        'O=C(NCC1CCCO1)Nc1ccccc1',  # Correct aniline
        'O=C(NCC1CCCO1)NC1CCCCC1',  # Wrong cyclohexyl
    ]

    for smi in test_smiles:
        try:
            # Get Molcart molfile from SMILES
            cur.execute(f"SELECT CONVERT(bin2mol(mol2bin('{smi}', 'smiles')) USING UTF8)")
            molcart_molfile = cur.fetchall()[0][0]
            # Get RDKit SMILES from that Molcart molfile
            if molcart_molfile:
                rdkit_smi = rdkit_standardize(molcart_molfile)
            else:
                rdkit_smi = "NULL"
            changed = "CHANGED!" if smi != rdkit_smi else "OK"
            print(f"  {smi:45s} -> Molcart molfile -> RDKit: {rdkit_smi:45s} [{changed}]")
        except Exception as e:
            print(f"  {smi}: ERROR: {e}")


    print("\n\n")
    print("=" * 120)
    print("TEST 6: Timeline - check all batches for affected compounds, with their")
    print("  chem_info compound_id and chem_info molfile SMILES")
    print("=" * 120)

    for cid in affected_compounds:
        print(f"\n  --- {cid} ---")
        cur.execute(f"""SELECT ci.regno, ci.jpage, ci.rdate, ci.compound_id, ci.sdfile_sequence
                       FROM chem_reg.chem_info ci
                       WHERE ci.compound_id = '{cid}'
                       ORDER BY ci.regno ASC""")
        rows = cur.fetchall()
        for r in rows:
            print(f"    regno={r[0]}, jpage={r[1]}, rdate={r[2]}, chem_info.compound_id={r[3]}, sdfile_seq={r[4]}")


    print("\n\n")
    print("=" * 120)
    print("TEST 7: Check if smiles_std_string differs from smiles_std")
    print("  (smiles_std_string was added in commit b8615b5, 2026-03-05)")
    print("=" * 120)

    for cid in affected_compounds:
        cur.execute(f"""SELECT smiles_std, smiles_std_string
                       FROM bcpvs.compound WHERE compound_id = '{cid}'""")
        rows = cur.fetchall()
        if rows:
            s1, s2 = rows[0]
            match = "SAME" if s1 == s2 else "DIFFERENT"
            print(f"  {cid}: smiles_std={s1}, smiles_std_string={s2}, {match}")


    print("\n\n")
    print("=" * 120)
    print("TEST 8: Check Molcart SMILES generation (bin2smiles) for the STORED smiles")
    print("  Does Molcart preserve aromatic/aliphatic distinction?")
    print("=" * 120)

    for smi in ['NC(=O)C1CCNCC1', 'NC(=O)c1ccncc1']:
        try:
            cur.execute(f"SELECT bin2smiles(mol2bin('{smi}', 'smiles'), 'mol')")
            result = cur.fetchall()[0][0]
            if result:
                result = result.decode() if isinstance(result, bytes) else result
            print(f"  bin2smiles(mol2bin('{smi}', 'smiles')) = {result}")
        except Exception as e:
            print(f"  {smi}: ERROR: {e}")


    cur.close()
    db.close()


if __name__ == '__main__':
    main()
