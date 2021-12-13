import mysql.connector
import mysql
import config

db_connection = mysql.connector.connect(
    host=config.fromOra['host'],
    user=config.fromOra['user'],
    passwd=config.fromOra['password']
)
cur = db_connection.cursor()


def dropTable(sTab):
    try:
        cur.execute("drop table " + sTab)
    except Exception as e:
        print(e)

allTabs = True
#allTabs = False

dropTable("hive.user_details")
dropTable("hive.project_details")


dropTable("chemspec.analytical_sum")
dropTable("chemspec.chemspec_salt")
dropTable("chemspec.chns_log")
dropTable("chemspec.color_tbl")
dropTable("chemspec.experiment")
dropTable("chemspec.ftir_log")
dropTable("chemspec.hplc_uv_log")
dropTable("chemspec.ms_log")
dropTable("chemspec.nmr_log")
dropTable("chemspec.other_log")
dropTable("chemspec.solvent_tbl")
dropTable("chemspec.spec_analysis")
dropTable("chemspec.storage_tbl")
dropTable("chemspec.chem_info")


dropTable("sfl_assay.rnai")

dropTable("microtube.matrix")
dropTable("microtube.matrix_position")
dropTable("microtube.matrix_tube")
dropTable("microtube.tube")
dropTable("microtube.tube_changes")
#dropTable("microtube.tube_log")

dropTable("loctree.location_access")
dropTable("loctree.location_hierarchy")
dropTable("loctree.location_type")
dropTable("loctree.locations")
dropTable("loctree.subpositions")


dropTable("cool.config")
dropTable("cool.config_key")
dropTable("cool.map96to384")
dropTable("cool.plate")
dropTable("cool.plate_type")
dropTable("cool.plating_sequence")
dropTable("cool.solubility_problem")



dropTable("assay.cerep_assays")
dropTable("assay.cerep_functional_dr")
dropTable("assay.cerep_functional_sp")
dropTable("assay.cerep_ki")
dropTable("assay.cerep_receptors")
dropTable("assay.cerep_screen")
dropTable("assay.detection_types")
dropTable("assay.lcb_dr")
dropTable("assay.lcb_dr_invalid")
dropTable("assay.lcb_sp_schneider_cysm")
dropTable("assay.lcb_sp_schneider_cysm_invalid")
dropTable("assay.metabolite")
dropTable("assay.upstate_sp")

dropTable("assay.ba_mds_dr")
dropTable("assay.ba_biofocus_sp")
dropTable("assay.ba_biofocus_dr")
dropTable("assay.ba_bf_camp_sp")
dropTable("assay.ba_bf_camp_dr")
dropTable("assay.assay_types")
dropTable("assay.lcb_sp")
dropTable("assay.lcb_sp_invalid")

dropTable("glass.globals")
dropTable("glass.vial_log")
dropTable("glass.vial")
dropTable("glass.hive_stats")


dropTable("screen.aa_cyp3a4_stab")
dropTable("screen.aa_elph_1")
dropTable("screen.aa_ez4u_screen")
dropTable("screen.aa_induction")
dropTable("screen.aa_ips")
dropTable("screen.aa_ltssolubility")
dropTable("screen.aa_metstab")
dropTable("screen.aa_mmt_ames")
dropTable("screen.aa_p450inhib")
dropTable("screen.aa_permeabilitet")
dropTable("screen.aa_pka")
dropTable("screen.aa_protein_binding")
dropTable("screen.aa_solubility")
dropTable("screen.aa_spec_analysis")
dropTable("screen.ba_binding_sp")
dropTable("screen.ba_ca_ic50")
dropTable("screen.ba_ca_sp")
dropTable("screen.ba_camp_efficacy")
dropTable("screen.ba_cell_sp")
dropTable("screen.ba_enz_inh_comp")
dropTable("screen.ba_enzyme_sp")
dropTable("screen.ba_functional_dr")
dropTable("screen.ba_herg_dr")
dropTable("screen.ba_herg_sp")
dropTable("screen.ba_hit_88")
dropTable("screen.ba_hts_80_2")
dropTable("screen.ba_radio_sp")
dropTable("screen.ba_rec_ec50")
dropTable("screen.hts_testsets")


dropTable("bcpvs.jcmol_moltable")
dropTable("bcpvs.compound_library")
dropTable("bcpvs.compound_suppliers")
dropTable("bcpvs.compound_type")
dropTable("bcpvs.cryst_solvent")
dropTable("bcpvs.external_id")
dropTable("bcpvs.external_id_invalid")
dropTable("bcpvs.id_ok")
dropTable("bcpvs.ip_rights_values")
dropTable("bcpvs.lcb_batch")
dropTable("bcpvs.molmass")
dropTable("bcpvs.periodic_table")
dropTable("bcpvs.product_type")
dropTable("bcpvs.screening_suppliers")
dropTable("bcpvs.structure_change_log")


if allTabs:
    dropTable("bcpvs.batch")
    dropTable("bcpvs.batch_invalid")
    dropTable("bcpvs.compound")
