#!/usr/bin/env python3
"""
update_structures.py - Dry run mode

Reads chemreg_tofix_claude.csv and analyzes what updates would be needed
to make bcpvs.compound match the supplier structure from chem_reg.chem_info.

For each compound_id:
  1. Gets current bcpvs.compound structure data
  2. Finds ALL batches pointing to that compound_id
  3. For each batch, fetches the MOLFILE from chem_reg.chem_info or chemspec.chem_info
  4. Standardizes each MOLFILE to canonical SMILES
  5. Checks if all batches agree on the supplier structure
  6. Categorizes as SIMPLE_UPDATE (safe) or CONFLICT (needs manual review)

Outputs a detailed report CSV: proposed_updates.csv
"""

import csv
import sys
import MySQLdb
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.warning')
RDLogger.DisableLog('rdApp.error')

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

import config

# ---------------------------------------------------------------------------
# Chemistry helpers
# ---------------------------------------------------------------------------

def standardize_smiles(smiles):
    """Standardize a SMILES: cleanup, fragment parent, uncharge, canonicalize tautomer."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    try:
        params = rdMolStandardize.CleanupParameters()
        params.tautomerRemoveSp3Stereo = False
        params.tautomerRemoveBondStereo = False
        params.tautomerRemoveIsotopicHs = False
        clean_mol = rdMolStandardize.Cleanup(mol, params)
        parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
        uncharged_mol = rdMolStandardize.Uncharger().uncharge(parent_mol)
        te = rdMolStandardize.TautomerEnumerator(params)
        canon_mol = te.Canonicalize(uncharged_mol)
        return Chem.MolToSmiles(canon_mol), canon_mol
    except Exception as e:
        logger.warning(f"Standardization failed for {smiles}: {e}")
        return None, None


def standardize_molfile(molfile):
    """Standardize a molfile: cleanup, fragment parent, uncharge, canonicalize tautomer."""
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    if mol is None:
        return None, None
    try:
        params = rdMolStandardize.CleanupParameters()
        params.tautomerRemoveSp3Stereo = False
        params.tautomerRemoveBondStereo = False
        params.tautomerRemoveIsotopicHs = False
        clean_mol = rdMolStandardize.Cleanup(mol, params)
        parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
        uncharged_mol = rdMolStandardize.Uncharger().uncharge(parent_mol)
        te = rdMolStandardize.TautomerEnumerator(params)
        canon_mol = te.Canonicalize(uncharged_mol)
        return Chem.MolToSmiles(canon_mol), canon_mol
    except Exception as e:
        logger.warning(f"Standardization of molfile failed: {e}")
        return None, None


def compute_properties(mol):
    """Compute molecular properties from an RDKit mol object."""
    if mol is None:
        return {}
    smiles = Chem.MolToSmiles(mol)
    inchi = MolToInchi(mol)
    inchi_key = InchiToInchiKey(inchi) if inchi else None
    mf = rdMolDescriptors.CalcMolFormula(mol)
    monoiso_mass = Descriptors.ExactMolWt(mol)
    return {
        'smiles': smiles,
        'inchi': inchi,
        'inchi_key': inchi_key,
        'mf': mf,
        'monoiso_mass': round(monoiso_mass, 6),
    }


# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------

def get_db_connection():
    return MySQLdb.connect(
        host=config.database['host'],
        user=config.database['user'],
        passwd=config.database['password'],
        db=config.database['db'],
        charset='utf8mb4'
    )


def get_compound_info(cur, compound_id):
    """Get current structure info from bcpvs.compound."""
    cur.execute("""
        SELECT compound_id, smiles_std, smiles_std_string, mf,
               SEP_MOL_MONOISO_MASS, inchi, inchi_key
        FROM bcpvs.compound
        WHERE compound_id = %s
    """, (compound_id,))
    row = cur.fetchone()
    if row:
        return {
            'compound_id': row[0],
            'smiles_std': row[1],
            'smiles_std_string': row[2],
            'mf': row[3],
            'monoiso_mass': row[4],
            'inchi': row[5],
            'inchi_key': row[6],
        }
    return None


def get_all_batches(cur, compound_id):
    """Get all batches pointing to a compound_id."""
    cur.execute("""
        SELECT notebook_ref, compound_id, chemreg_regno, chemspec_regno
        FROM bcpvs.batch
        WHERE compound_id = %s
    """, (compound_id,))
    rows = cur.fetchall()
    return [
        {
            'notebook_ref': r[0],
            'compound_id': r[1],
            'chemreg_regno': r[2],
            'chemspec_regno': r[3],
        }
        for r in rows
    ]


def get_molfile_from_chemreg(cur, regno):
    """Get MOLFILE from chem_reg.chem_info."""
    cur.execute("""
        SELECT MOLFILE FROM chem_reg.chem_info WHERE regno = %s
    """, (regno,))
    row = cur.fetchone()
    if row and row[0]:
        return row[0]
    return None


def get_molfile_from_chemspec(cur, regno):
    """Get molfile from chemspec.chem_info."""
    cur.execute("""
        SELECT molfile FROM chemspec.chem_info WHERE regno = %s
    """, (regno,))
    row = cur.fetchone()
    if row and row[0]:
        return row[0]
    return None


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def analyze_compound(cur, compound_id, csv_rows):
    """
    Analyze a compound_id and determine what update action is needed.

    Returns a dict with analysis results.
    """
    result = {
        'compound_id': compound_id,
        'action': None,           # SIMPLE_UPDATE, CONFLICT, ERROR, SKIP
        'reason': '',
        'current': None,          # current bcpvs.compound properties
        'proposed': None,         # proposed new properties
        'total_batches': 0,
        'batch_details': [],      # per-batch info
        'csv_batches': len(csv_rows),
    }

    # Get current compound info
    compound_info = get_compound_info(cur, compound_id)
    if not compound_info:
        result['action'] = 'ERROR'
        result['reason'] = f'compound_id {compound_id} not found in bcpvs.compound'
        return result
    result['current'] = compound_info

    # Get all batches for this compound
    batches = get_all_batches(cur, compound_id)
    result['total_batches'] = len(batches)

    if len(batches) == 0:
        result['action'] = 'ERROR'
        result['reason'] = 'No batches found for compound_id'
        return result

    # For each batch, get the original supplier MOLFILE and standardize it
    batch_structures = {}  # notebook_ref -> standardized SMILES
    batch_details = []

    for batch in batches:
        nb_ref = batch['notebook_ref']
        detail = {
            'notebook_ref': nb_ref,
            'chemreg_regno': batch['chemreg_regno'],
            'chemspec_regno': batch['chemspec_regno'],
            'source': None,
            'molfile_found': False,
            'std_smiles': None,
        }

        molfile = None
        if batch['chemreg_regno']:
            detail['source'] = 'chem_reg'
            molfile = get_molfile_from_chemreg(cur, batch['chemreg_regno'])
        elif batch['chemspec_regno']:
            detail['source'] = 'chemspec'
            molfile = get_molfile_from_chemspec(cur, batch['chemspec_regno'])
        else:
            detail['source'] = 'NONE'
            detail['molfile_found'] = False
            batch_details.append(detail)
            continue

        if molfile:
            detail['molfile_found'] = True
            std_smiles, std_mol = standardize_molfile(molfile)
            detail['std_smiles'] = std_smiles
            if std_smiles:
                batch_structures[nb_ref] = std_smiles
        else:
            detail['molfile_found'] = False

        batch_details.append(detail)

    result['batch_details'] = batch_details

    # Check if we could get structures for all batches
    batches_with_structures = [d for d in batch_details if d['std_smiles']]
    batches_without_structures = [d for d in batch_details if d['molfile_found'] and not d['std_smiles']]
    batches_without_molfile = [d for d in batch_details if not d['molfile_found']]

    if not batches_with_structures:
        result['action'] = 'ERROR'
        result['reason'] = 'Could not standardize any batch MOLFILE'
        return result

    # Get unique standardized SMILES across all batches
    unique_smiles = set(batch_structures.values())

    if len(unique_smiles) == 1:
        # All batches agree on the structure - simple update
        target_smiles = unique_smiles.pop()
        target_mol = Chem.MolFromSmiles(target_smiles)
        if target_mol is None:
            result['action'] = 'ERROR'
            result['reason'] = f'Could not parse agreed-upon SMILES: {target_smiles}'
            return result

        proposed = compute_properties(target_mol)

        # Check if compound already has the correct structure
        current_std, _ = standardize_smiles(compound_info['smiles_std']) if compound_info['smiles_std'] else (None, None)
        if current_std == target_smiles:
            result['action'] = 'SKIP'
            result['reason'] = 'bcpvs.compound already has the correct standardized structure'
            result['proposed'] = proposed
            return result

        result['action'] = 'SIMPLE_UPDATE'
        result['reason'] = (
            f'{len(batches)} batch(es), all agree on supplier structure'
            + (f' ({len(batches_without_molfile)} batch(es) had no MOLFILE)' if batches_without_molfile else '')
        )
        result['proposed'] = proposed

    else:
        # Multiple distinct structures across batches - conflict
        result['action'] = 'CONFLICT'
        smiles_list = [f"{d['notebook_ref']}={d['std_smiles']}" for d in batch_details if d['std_smiles']]
        result['reason'] = (
            f'{len(batches)} batches with {len(unique_smiles)} distinct structures: '
            + '; '.join(smiles_list)
        )
        # Still compute what the CSV rows suggest as the correct structure
        # Use the supplier_smiles from the first CSV row
        supplier_smiles = csv_rows[0].get('supplier_smiles', '')
        if supplier_smiles:
            std_supplier, std_mol = standardize_smiles(supplier_smiles)
            if std_mol:
                result['proposed'] = compute_properties(std_mol)

    return result


def main():
    # Read CSV file
    csv_file = 'chemreg_tofix_claude.csv'
    logger.info(f'Reading {csv_file}...')

    rows_by_compound = defaultdict(list)
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows_by_compound[row['compound_id']].append(row)

    logger.info(f'Found {sum(len(v) for v in rows_by_compound.values())} rows for '
                f'{len(rows_by_compound)} unique compound_ids')

    # Connect to database
    conn = get_db_connection()
    cur = conn.cursor()

    # Analyze each compound
    results = []
    counts = defaultdict(int)

    for i, (compound_id, csv_rows) in enumerate(sorted(rows_by_compound.items()), 1):
        if i % 50 == 0:
            logger.info(f'Processing {i}/{len(rows_by_compound)}...')

        analysis = analyze_compound(cur, compound_id, csv_rows)
        results.append(analysis)
        counts[analysis['action']] += 1

    cur.close()
    conn.close()

    # Print summary
    logger.info('=' * 60)
    logger.info('ANALYSIS SUMMARY')
    logger.info('=' * 60)
    logger.info(f"Total compound_ids analyzed: {len(results)}")
    for action in ['SIMPLE_UPDATE', 'CONFLICT', 'SKIP', 'ERROR']:
        logger.info(f"  {action}: {counts.get(action, 0)}")

    # Write detailed report
    report_file = 'proposed_updates.csv'
    with open(report_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'compound_id',
            'action',
            'reason',
            'total_batches',
            'csv_batches',
            'current_smiles_std',
            'current_mf',
            'current_monoiso_mass',
            'current_inchi_key',
            'proposed_smiles',
            'proposed_mf',
            'proposed_monoiso_mass',
            'proposed_inchi',
            'proposed_inchi_key',
            'batch_details',
        ])

        for r in results:
            current = r['current'] or {}
            proposed = r['proposed'] or {}

            # Format batch details as compact string
            bd_parts = []
            for bd in r['batch_details']:
                source = bd.get('source', '?')
                regno = bd.get('chemreg_regno') or bd.get('chemspec_regno') or '?'
                smiles = bd.get('std_smiles', 'N/A')
                mf = 'yes' if bd.get('molfile_found') else 'no'
                bd_parts.append(f"{bd['notebook_ref']}[{source}:{regno},mf={mf},smi={smiles}]")
            batch_str = ' | '.join(bd_parts)

            writer.writerow([
                r['compound_id'],
                r['action'],
                r['reason'],
                r['total_batches'],
                r['csv_batches'],
                current.get('smiles_std', ''),
                current.get('mf', ''),
                current.get('monoiso_mass', ''),
                current.get('inchi_key', ''),
                proposed.get('smiles', ''),
                proposed.get('mf', ''),
                proposed.get('monoiso_mass', ''),
                proposed.get('inchi', ''),
                proposed.get('inchi_key', ''),
                batch_str,
            ])

    logger.info(f'Report written to {report_file}')

    # Also write conflict details to a separate file for easier review
    conflict_file = 'conflicts.csv'
    conflicts = [r for r in results if r['action'] == 'CONFLICT']
    if conflicts:
        with open(conflict_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'compound_id',
                'total_batches',
                'notebook_ref',
                'source',
                'regno',
                'std_smiles',
                'molfile_found',
            ])
            for r in conflicts:
                for bd in r['batch_details']:
                    writer.writerow([
                        r['compound_id'],
                        r['total_batches'],
                        bd['notebook_ref'],
                        bd.get('source', ''),
                        bd.get('chemreg_regno') or bd.get('chemspec_regno') or '',
                        bd.get('std_smiles', ''),
                        bd.get('molfile_found', ''),
                    ])
        logger.info(f'Conflict details written to {conflict_file}')


if __name__ == '__main__':
    main()
