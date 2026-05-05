#!/usr/bin/env python3
"""Backfill the inchi_key column on bcpvs[_test].compound and chem_reg[_test].chem_info.

Run after applying 2026-04-29-add-chem_info-inchi_key.sql.

For bcpvs.compound rows where inchi_key IS NULL, parses smiles_std and
computes the InChIKey. For chem_reg.chem_info rows where inchi_key IS NULL
and a molfile is present, parses the molfile and computes the InChIKey
*after* the same standardization pipeline as the live registration code
(cleanup -> fragment parent -> uncharge; tautomer canonicalization is
deliberately skipped because it is the source of the non-determinism we
are trying to remove from the registration path).

Usage:
    python3 backfill_inchi_key.py --db bcpvs_test           # dry run
    python3 backfill_inchi_key.py --db bcpvs_test --apply
    python3 backfill_inchi_key.py --db chem_reg_test --apply
    python3 backfill_inchi_key.py --db bcpvs --apply
    python3 backfill_inchi_key.py --db chem_reg --apply
"""
import argparse
import logging
import sys

import MySQLdb
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi

RDLogger.DisableLog('rdApp.*')

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def _standardize_no_tautomer(mol):
    """Cleanup + FragmentParent + Uncharger. No tautomer canonicalization
    (InChI handles tautomer normalization itself, and the RDKit tautomer
    enumerator is the source of the non-determinism we are working
    around)."""
    if mol is None:
        return None
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    clean_mol = rdMolStandardize.Cleanup(mol, params)
    parent_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    return rdMolStandardize.Uncharger().uncharge(parent_mol)


def _inchi_key_from_mol(mol):
    if mol is None:
        return None
    inchi = MolToInchi(mol, options='-w')
    if not inchi:
        return None
    key = InchiToInchiKey(inchi)
    return key or None


def inchi_key_from_smiles(smiles):
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        return _inchi_key_from_mol(_standardize_no_tautomer(mol))
    except Exception as e:
        logger.warning(f"smiles failed: {smiles[:80]}... {e}")
        return None


def inchi_key_from_molfile(molfile):
    if not molfile:
        return None
    try:
        mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
        return _inchi_key_from_mol(_standardize_no_tautomer(mol))
    except Exception as e:
        logger.warning(f"molfile failed: {e}")
        return None


def backfill_compound(cur, db, apply_changes):
    cur.execute(f"""SELECT compound_id, smiles_std FROM {db}.compound
                    WHERE inchi_key IS NULL OR inchi_key = ''""")
    rows = cur.fetchall()
    logger.info(f"{db}.compound: {len(rows)} rows missing inchi_key")
    updated = 0
    failed = 0
    for compound_id, smiles in rows:
        key = inchi_key_from_smiles(smiles)
        if not key:
            failed += 1
            continue
        if apply_changes:
            cur.execute(f"UPDATE {db}.compound SET inchi_key = %s WHERE compound_id = %s",
                        (key, compound_id))
        updated += 1
    logger.info(f"{db}.compound: {updated} backfilled, {failed} failed")


def backfill_chem_info(cur, db, apply_changes):
    cur.execute(f"""SELECT regno, MOLFILE FROM {db}.chem_info
                    WHERE (inchi_key IS NULL OR inchi_key = '')
                      AND MOLFILE IS NOT NULL""")
    rows = cur.fetchall()
    logger.info(f"{db}.chem_info: {len(rows)} rows missing inchi_key")
    updated = 0
    failed = 0
    for regno, molfile in rows:
        key = inchi_key_from_molfile(molfile)
        if not key:
            failed += 1
            continue
        if apply_changes:
            cur.execute(f"UPDATE {db}.chem_info SET inchi_key = %s WHERE regno = %s",
                        (key, regno))
        updated += 1
    logger.info(f"{db}.chem_info: {updated} backfilled, {failed} failed")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--db', required=True,
                   choices=['bcpvs', 'bcpvs_test', 'chem_reg', 'chem_reg_test'])
    p.add_argument('--apply', action='store_true', help='Without this flag run as a dry run.')
    p.add_argument('--config', default=None, help='Path to a config.py with database dict.')
    args = p.parse_args()

    if args.config:
        sys.path.insert(0, args.config)
    import config

    conn = MySQLdb.connect(host=config.database['host'],
                           user=config.database['user'],
                           passwd=config.database['password'],
                           charset='utf8mb4', use_unicode=True)
    conn.autocommit(True)
    cur = conn.cursor()

    if args.db.startswith('bcpvs'):
        backfill_compound(cur, args.db, args.apply)
    else:
        backfill_chem_info(cur, args.db, args.apply)


if __name__ == '__main__':
    main()
