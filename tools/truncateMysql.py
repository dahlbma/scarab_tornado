import mysql.connector
import mysql
import config

db_connection = mysql.connector.connect(
    host=config.fromOra['host'],
    user=config.fromOra['user'],
    passwd=config.fromOra['password']
)
cur = db_connection.cursor()
cur.execute('SET FOREIGN_KEY_CHECKS = 0')

def truncateTable(sTab):
    try:
        cur.execute("truncate table " + sTab)
    except Exception as e:
        print(e)

allTabs = True
#allTabs = False

truncateTable("chem_reg.CHEM")
truncateTable("chem_reg.CHEM_MOL")
truncateTable("chem_reg.CHEM_MOL_key")
truncateTable("chem_reg.CHEM_MOL_keysim")
truncateTable("chem_reg.chem_info")
truncateTable("chem_reg.chem_info_mol_key")
truncateTable("chem_reg.chem_info_mol_ukey")
truncateTable("chem_reg.chemreg_dist")
truncateTable("chem_reg.tmp_mol")
<<<<<<< HEAD

=======
>>>>>>> 5ffcb722de566b5de2b190e044f06ddec9ffe466
quit()

truncateTable("bcpvs.JCMOL_MOLTABLE_MOL_keysim")
truncateTable("bcpvs.JCMOL_MOLTABLE_MOL_key")
truncateTable("bcpvs.JCMOL_MOLTABLE_MOL")
truncateTable("bcpvs.JCMOL_MOLTABLE")
truncateTable("bcpvs.compound_library")
truncateTable("bcpvs.compound_suppliers")
truncateTable("bcpvs.compound_type")
truncateTable("bcpvs.cryst_solvent")
truncateTable("bcpvs.external_id")
truncateTable("bcpvs.external_id_invalid")
truncateTable("bcpvs.id_ok")
truncateTable("bcpvs.ip_rights_values")
truncateTable("bcpvs.lcb_batch")
truncateTable("bcpvs.molmass")
truncateTable("bcpvs.periodic_table")
truncateTable("bcpvs.product_type")
truncateTable("bcpvs.screening_suppliers")
truncateTable("bcpvs.structure_change_log")


truncateTable("chemspec.analytical_sum")
truncateTable("chemspec.chemspec_salt")
truncateTable("chemspec.chns_log")
truncateTable("chemspec.color_tbl")
truncateTable("chemspec.experiment")
truncateTable("chemspec.ftir_log")
truncateTable("chemspec.hplc_uv_log")
truncateTable("chemspec.ms_log")
truncateTable("chemspec.nmr_log")
truncateTable("chemspec.other_log")
truncateTable("chemspec.solvent_tbl")
truncateTable("chemspec.spec_analysis")
truncateTable("chemspec.storage_tbl")
truncateTable("chemspec.chem_info")


truncateTable("sfl_assay.rnai")

truncateTable("microtube.matrix_tube_invalid")
truncateTable("microtube.matrix_tube")
truncateTable("microtube.matrix")
truncateTable("microtube.matrix_position")
truncateTable("microtube.tube")
truncateTable("microtube.tube_invalid")
truncateTable("microtube.tube_changes")
#truncateTable("microtube.tube_log")

truncateTable("loctree.location_access")
truncateTable("loctree.location_hierarchy")
truncateTable("loctree.location_type")
truncateTable("loctree.locations")
truncateTable("loctree.subpositions")


truncateTable("cool.config")
truncateTable("cool.config_key")
truncateTable("cool.map96to384")
truncateTable("cool.plate")
truncateTable("cool.plate_type")
truncateTable("cool.plating_sequence")
truncateTable("cool.solubility_problem")


truncateTable("assay.cerep_assays")
truncateTable("assay.cerep_functional_dr")
truncateTable("assay.cerep_functional_sp")
truncateTable("assay.cerep_ki")
truncateTable("assay.cerep_receptors")
truncateTable("assay.cerep_screen")
truncateTable("assay.detection_types")
truncateTable("assay.lcb_dr")
truncateTable("assay.lcb_dr_invalid")
truncateTable("assay.lcb_sp_schneider_cysm")
truncateTable("assay.lcb_sp_schneider_cysm_invalid")
truncateTable("assay.metabolite")
truncateTable("assay.upstate_sp")

truncateTable("assay.ba_mds_dr")
truncateTable("assay.ba_biofocus_sp")
truncateTable("assay.ba_biofocus_dr")
truncateTable("assay.ba_bf_camp_sp")
truncateTable("assay.ba_bf_camp_dr")
truncateTable("assay.assay_types")
truncateTable("assay.lcb_sp")
truncateTable("assay.lcb_sp_invalid")

truncateTable("glass.globals")
truncateTable("glass.vial_log")
truncateTable("glass.vial_log_invalid")
truncateTable("glass.vial")
truncateTable("glass.vial_invalid")
truncateTable("glass.hive_stats")


truncateTable("screen.aa_cyp3a4_stab")
truncateTable("screen.aa_elph_1")
truncateTable("screen.aa_ez4u_screen")
truncateTable("screen.aa_induction")
truncateTable("screen.aa_ips")
truncateTable("screen.aa_ltssolubility")
truncateTable("screen.aa_metstab")
truncateTable("screen.aa_mmt_ames")
truncateTable("screen.aa_p450inhib")
truncateTable("screen.aa_permeabilitet")
truncateTable("screen.aa_pka")
truncateTable("screen.aa_protein_binding")
truncateTable("screen.aa_solubility")
truncateTable("screen.aa_spec_analysis")
truncateTable("screen.ba_binding_sp")
truncateTable("screen.ba_ca_ic50")
truncateTable("screen.ba_ca_sp")
truncateTable("screen.ba_camp_efficacy")
truncateTable("screen.ba_cell_sp")
truncateTable("screen.ba_enz_inh_comp")
truncateTable("screen.ba_enzyme_sp")
truncateTable("screen.ba_functional_dr")
truncateTable("screen.ba_herg_dr")
truncateTable("screen.ba_herg_sp")
truncateTable("screen.ba_hit_88")
truncateTable("screen.ba_hts_80_2")
truncateTable("screen.ba_radio_sp")
truncateTable("screen.ba_rec_ec50")
truncateTable("screen.hts_testsets")


if allTabs:
    truncateTable("bcpvs.batch")
    truncateTable("bcpvs.batch_invalid")
    truncateTable("bcpvs.compound")


cur.execute('SET FOREIGN_KEY_CHECKS = 1')
