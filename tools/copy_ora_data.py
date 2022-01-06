import cx_Oracle
import mysql.connector
import mysql
import pandas
import pandas_oracle.tools as pt
import sqlalchemy
import config

con = cx_Oracle.connect(config.Ora['user'],
                        config.Ora['password'],
                        config.Ora['host'])
print("Database version:", con.version)

db_connection = mysql.connector.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password']
)

cur = db_connection.cursor()
## opening conn
#conn = pt.open_connection("config.yml")

#query1 = "select * from  hive.hive"
## passing the conn object to the query_to_df
#df1 = pt.query_to_df(query1, conn, 10000)

#engineBEEHIVE = sqlalchemy.create_engine('mysql+pymysql://username:password@localhost/SCHEMA')
engineHIVE = sqlalchemy.create_engine('mysql+pymysql://' + config.HIVE)
engineBCPVS = sqlalchemy.create_engine('mysql+pymysql://' + config.BCPVS)
engineASSAY = sqlalchemy.create_engine('mysql+pymysql://' + config.ASSAY)
engineCOOL = sqlalchemy.create_engine('mysql+pymysql://' + config.COOL)
engineGLASS = sqlalchemy.create_engine('mysql+pymysql://' + config.GLASS)
engineSFL_ASSAY = sqlalchemy.create_engine('mysql+pymysql://' + config.SFL_ASSAY)
engineLOCTREE = sqlalchemy.create_engine('mysql+pymysql://' + config.LOCTREE)
engineMICROTUBE = sqlalchemy.create_engine('mysql+pymysql://' + config.MICROTUBE)
engineCHEMSPEC = sqlalchemy.create_engine('mysql+pymysql://' + config.CHEMSPEC)
engineSCREEN = sqlalchemy.create_engine('mysql+pymysql://' + config.SCREEN)
engineSCREENMOL = sqlalchemy.create_engine('mysql+pymysql://' + config.SCREENMOL)

################################################

def copyTable(sEngine, sTable, sName):
    tableExists = False
    query = "select * from " + sTable
    try:
        cur.execute(query)
        myresult = cur.fetchall()
        tableExists = True
    except:
        tableExists = False

    if tableExists:
        return False
    
    query = "select * from " + sTable
    df2 = pandas.read_sql(sql=query, 
                          con=con)
    bReturn = False
    try:
        print('copy ' + sTable)
        df2.to_sql(sName, sEngine, index=False,
        dtype={'MOLFILE':  sqlalchemy.types.Text()})
        sSchemaTab = sTable.split('.')[0] + '.' + sName
        #cur.execute(f"alter table {sSchemaTab} convert to character set latin1 collate latin1_swedish_ci")
        bReturn = True
    except Exception as e:
        print('Error on creating: ' + sTable + ': ' + str(e))
    return bReturn

'''
df2.to_sql(sName, sEngine, dtype={'index':  sqlalchemy.types.INTEGER()})

dtype={'name_of_datefld': sqlalchemy.DateTime(), 
       'name_of_intfld': sqlalchemy.types.INTEGER(),
       'name_of_strfld': sqlalchemy.types.VARCHAR(length=30),
       'name_of_floatfld': sqlalchemy.types.Float(precision=3, asdecimal=True),
       'name_of_booleanfld': sqlalchemy.types.Boolean}
'''


################################################
# HIVE tables

if copyTable(engineHIVE, 'hive.user_details', 'user_details'):
    cur.execute("ALTER TABLE hive.user_details Modify column userid varchar(40)")
    cur.execute("ALTER TABLE hive.user_details Modify column firstname varchar(40)")
    cur.execute("ALTER TABLE hive.user_details Modify column lastname varchar(40)")
    cur.execute("ALTER TABLE hive.user_details Modify column email varchar(80)")
    cur.execute("ALTER TABLE hive.user_details Modify column unix_account varchar(4)")

    cur.execute("ALTER TABLE hive.user_details Modify column location varchar(30)")
    cur.execute("ALTER TABLE hive.user_details Modify column organisation varchar(50)")
    cur.execute("ALTER TABLE hive.user_details Modify column fullname varchar(82)")
    cur.execute("ALTER TABLE hive.user_details Modify column organization varchar(255)")
    cur.execute("""update hive.user_details set organization = 'chemistry' where userid in
                   ('ANGUST', 'MAHARA', 'SLAS', 'DAHLBMA', 'BISJO')""")
    
    try:
        cur.execute("""ALTER TABLE hive.user_details CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on hive.user_details.pkey')
##
if copyTable(engineHIVE, 'hive.project_details_lcb', 'project_details'):
    cur.execute("ALTER TABLE hive.project_details Modify column project_name varchar(40)")
    cur.execute("ALTER TABLE hive.project_details Modify column project_leader varchar(10)")

    try:
        cur.execute("""ALTER TABLE hive.project_details
                       CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on hive.project_details.pkey')

################################################
# BCPVS tables
if copyTable(engineBCPVS, 'bcpvs.compound', 'compound'):
    cur.execute("ALTER TABLE bcpvs.compound Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE bcpvs.compound Modify column suffix varchar(5)")
    cur.execute("ALTER TABLE bcpvs.compound Modify column ip_rights varchar(30)")
    cur.execute("ALTER TABLE bcpvs.compound Modify column mf varchar(200)")
    try:
        cur.execute("""CREATE UNIQUE INDEX compound_idx ON bcpvs.compound(compound_id)""")
    except Exception as e:
        print('Faile to create index on bcpvs.compound ' + str(e))

        cur.execute(f'''CREATE TABLE bcpvs.compound_id_sequence (
        pkey int NOT NULL,
        PRIMARY KEY (pkey))''')
        cur.execute(f'''insert into bcpvs.compound_id_sequence set pkey = 600000''')
 
    
##

if copyTable(engineBCPVS, 'bcpvs.batch', 'batch'):
    cur.execute("ALTER TABLE bcpvs.batch Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column notebook_ref varchar(16)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column submitter varchar(30)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column project varchar(30)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column supplier varchar(100)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column colour varchar(20)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column physical_state varchar(20)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column cryst_solvent varchar(20)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column library_id varchar(100)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column compound_type varchar(20)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column product_type varchar(20)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column supplier_id varchar(200)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column supplier_batch varchar(100)")
    cur.execute("ALTER TABLE bcpvs.batch Modify column chemspec_regno varchar(16)")

    cur.execute("ALTER TABLE bcpvs.batch ADD COLUMN suffix varchar(20)")
    
    """ These batches have 2 entries
    BC9614001
    BJ1804006
    BJ1804007
    BJ1804011
    BH1101002

    compound_id CBK27777 in batch doesn't exist in the compound table
    """
    cur.execute("""create table bcpvs.batch_invalid as 
    (select * from bcpvs.batch where notebook_ref in (
    'BC9614001','BJ1804006', 'BJ1804007', 'BJ1804011', 'BH1101002'
    ) or compound_id not in (select compound_id from bcpvs.compound))""")
    
    cur.execute("""delete from bcpvs.batch where notebook_ref in 
    ('BC9614001','BJ1804006', 'BJ1804007', 'BJ1804011', 'BH1101002')
    """)

    cur.execute("""delete from bcpvs.batch where compound_id in 
    (select * from
    (select distinct(compound_id) from bcpvs.batch where compound_id not in 
    (select compound_id from bcpvs.compound))tmpTbl)
    """)


    try:
        cur.execute("""CREATE UNIQUE INDEX batch_idx ON bcpvs.batch(notebook_ref)""")
    except Exception as e:
        print(str(e))
        print('Error creating index')

    cur.execute("""ALTER TABLE bcpvs.batch
                  ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

##
#if copyTable(engineBCPVS, 'bcpvs.jcmol_moltable', 'jcmol_moltable'):
#    cur.execute("ALTER TABLE bcpvs.jcmol_moltable Modify column compound_id varchar(26)")
#
#    try:
#        cur.execute("""ALTER TABLE bcpvs.jcmol_moltable CHANGE cd_id cd_id bigint AUTO_INCREMENT PRIMARY KEY""")
#    except Exception as e:
#        print('Faile to create index on bcpvs.jcmol_moltable ' + str(e))
#
#    cur.execute("""ALTER TABLE bcpvs.jcmol_moltable
#                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

##
if copyTable(engineBCPVS, 'bcpvs.compound_library', 'compound_library'):
    cur.execute("ALTER TABLE bcpvs.compound_library Modify column library_name varchar(50)")

    try:
        cur.execute("""ALTER TABLE bcpvs.compound_library CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print('Faile to create index on bcpvs.compound_library ' + str(e))
##
if copyTable(engineBCPVS, 'bcpvs.compound_suppliers', 'compound_suppliers'):
    cur.execute("ALTER TABLE bcpvs.compound_suppliers Modify column name varchar(40)")
##
if copyTable(engineBCPVS, 'bcpvs.compound_suppliers', 'compound_suppliers'):
    cur.execute("ALTER TABLE bcpvs.compound_suppliers Modify column name varchar(50)")
##
if copyTable(engineBCPVS, 'bcpvs.compound_type', 'compound_type'):
    cur.execute("ALTER TABLE bcpvs.compound_type Modify column type varchar(20)")
##
if copyTable(engineBCPVS, 'bcpvs.cryst_solvent', 'cryst_solvent'):
    cur.execute("ALTER TABLE bcpvs.cryst_solvent Modify column solvent varchar(20)")
##
if copyTable(engineBCPVS, 'bcpvs.external_id', 'external_id'):
    cur.execute("ALTER TABLE bcpvs.external_id Modify column external_id varchar(100)")
    cur.execute("ALTER TABLE bcpvs.external_id Modify column compound_id varchar(26)")

    cur.execute("""create table bcpvs.external_id_invalid as 
    (select * from bcpvs.external_id where compound_id not in (
    select compound_id from bcpvs.compound
    ))""")
    
    cur.execute("""delete from bcpvs.external_id where compound_id in 
    (select * from
    (select distinct(compound_id) from bcpvs.external_id where compound_id not in 
    (select compound_id from bcpvs.compound))tmpTbl)
    """)
    
    cur.execute("""ALTER TABLE bcpvs.external_id
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
    try:
        cur.execute("""ALTER TABLE bcpvs.external_id CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print('Faile to create index on bcpvs.external_id ' + str(e))
##
if copyTable(engineBCPVS, 'bcpvs.id_ok', 'id_ok'):
    cur.execute("ALTER TABLE bcpvs.id_ok Modify column ok varchar(10)")
##
if copyTable(engineBCPVS, 'bcpvs.ip_rights_values', 'ip_rights_values'):
    cur.execute("ALTER TABLE bcpvs.ip_rights_values Modify column ip_rights varchar(15)")
##
if copyTable(engineBCPVS, 'bcpvs.lcb_batch', 'lcb_batch'):
    cur.execute("ALTER TABLE bcpvs.lcb_batch Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE bcpvs.lcb_batch Modify column notebook_ref varchar(16)")
##
if copyTable(engineBCPVS, 'bcpvs.molmass', 'molmass'):
    cur.execute("ALTER TABLE bcpvs.molmass Modify column atom varchar(5)")
##
if copyTable(engineBCPVS, 'bcpvs.periodic_table', 'periodic_table'):
    cur.execute("ALTER TABLE bcpvs.periodic_table Modify column atom varchar(5)")
##
if copyTable(engineBCPVS, 'bcpvs.product_type', 'product_type'):
    cur.execute("ALTER TABLE bcpvs.product_type Modify column type varchar(30)")
##
if copyTable(engineBCPVS, 'bcpvs.screening_suppliers', 'screening_suppliers'):
    cur.execute("ALTER TABLE bcpvs.screening_suppliers Modify column name varchar(50)")

    try:
        cur.execute("""ALTER TABLE bcpvs.screening_suppliers CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print('Faile to create index on bcpvs.screening_suppliers ' + str(e))
##
if copyTable(engineBCPVS, 'bcpvs.structure_change_log', 'structure_change_log'):
    cur.execute("ALTER TABLE bcpvs.structure_change_log Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE bcpvs.structure_change_log
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

################################################
# SCREEN tables
if copyTable(engineSCREEN, 'screen.aa_cyp3a4_stab', 'aa_cyp3a4_stab'):
    cur.execute("ALTER TABLE screen.aa_cyp3a4_stab Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE screen.aa_cyp3a4_stab Modify column compound_batch varchar(20)")
    cur.execute("ALTER TABLE screen.aa_cyp3a4_stab Modify column project varchar(20)")

    cur.execute("""ALTER TABLE screen.aa_cyp3a4_stab
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_elph_1', 'aa_elph_1'):
    cur.execute("ALTER TABLE screen.aa_elph_1 Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE screen.aa_elph_1 Modify column compound_batch varchar(20)")
    cur.execute("ALTER TABLE screen.aa_elph_1 Modify column project varchar(20)")

    cur.execute("""ALTER TABLE screen.aa_elph_1
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_ez4u_screen', 'aa_ez4u_screen'):
    cur.execute("ALTER TABLE screen.aa_ez4u_screen Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_ez4u_screen
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_induction', 'aa_induction'):
    cur.execute("ALTER TABLE screen.aa_induction Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_induction
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_ips', 'aa_ips'):
    cur.execute("ALTER TABLE screen.aa_ips Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE screen.aa_ips Modify column compound_batch varchar(16)")

    cur.execute("""ALTER TABLE screen.aa_ips
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
##
if copyTable(engineSCREEN, 'screen.aa_ltssolubility', 'aa_ltssolubility'):
    cur.execute("ALTER TABLE screen.aa_ltssolubility Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_ltssolubility
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_metstab', 'aa_metstab'):
    cur.execute("ALTER TABLE screen.aa_metstab Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_metstab
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_mmt_ames', 'aa_mmt_ames'):
    cur.execute("ALTER TABLE screen.aa_mmt_ames Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_mmt_ames
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_p450inhib', 'aa_p450inhib'):
    cur.execute("ALTER TABLE screen.aa_p450inhib Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_p450inhib
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_permeabilitet', 'aa_permeabilitet'):
    cur.execute("ALTER TABLE screen.aa_permeabilitet Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_permeabilitet
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_pka', 'aa_pka'):
    cur.execute("ALTER TABLE screen.aa_pka Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_pka
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_protein_binding', 'aa_protein_binding'):
    cur.execute("ALTER TABLE screen.aa_protein_binding Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_protein_binding
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.aa_solubility', 'aa_solubility'):
    cur.execute("ALTER TABLE screen.aa_solubility Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE screen.aa_solubility Modify column compound_batch varchar(16)")

    cur.execute("""ALTER TABLE screen.aa_solubility
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
##
if copyTable(engineSCREEN, 'screen.aa_spec_analysis', 'aa_spec_analysis'):
    cur.execute("ALTER TABLE screen.aa_spec_analysis Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.aa_spec_analysis
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_binding_sp', 'ba_binding_sp'):
    cur.execute("ALTER TABLE screen.ba_binding_sp Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_binding_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_ca_ic50', 'ba_ca_ic50'):
    cur.execute("ALTER TABLE screen.ba_ca_ic50 Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_ca_ic50
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_ca_sp', 'ba_ca_sp'):
    cur.execute("ALTER TABLE screen.ba_ca_sp Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_ca_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_camp_efficacy', 'ba_camp_efficacy'):
    cur.execute("ALTER TABLE screen.ba_camp_efficacy Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_camp_efficacy
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_cell_sp', 'ba_cell_sp'):
    cur.execute("ALTER TABLE screen.ba_cell_sp Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_cell_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_enz_inh_comp', 'ba_enz_inh_comp'):
    cur.execute("ALTER TABLE screen.ba_enz_inh_comp Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_enz_inh_comp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_enzyme_sp', 'ba_enzyme_sp'):
    cur.execute("ALTER TABLE screen.ba_enzyme_sp Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_enzyme_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_functional_dr', 'ba_functional_dr'):
    cur.execute("ALTER TABLE screen.ba_functional_dr Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_functional_dr
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_herg_dr', 'ba_herg_dr'):
    cur.execute("ALTER TABLE screen.ba_herg_dr Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_herg_dr
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_herg_sp', 'ba_herg_sp'):
    cur.execute("ALTER TABLE screen.ba_herg_sp Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_herg_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_hit_88', 'ba_hit_88'):
    cur.execute("ALTER TABLE screen.ba_hit_88 Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_hit_88
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_hts_80_2', 'ba_hts_80_2'):
    cur.execute("ALTER TABLE screen.ba_hts_80_2 Modify column compound_id varchar(26)")

    cur.execute("""ALTER TABLE screen.ba_hts_80_2
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")
##
if copyTable(engineSCREEN, 'screen.ba_radio_sp', 'ba_radio_sp'):
    cur.execute("ALTER TABLE screen.ba_radio_sp Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE screen.ba_radio_sp Modify column compound_batch varchar(16)")

    cur.execute("""ALTER TABLE screen.ba_radio_sp
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
##
if copyTable(engineSCREEN, 'screen.ba_rec_ec50', 'ba_rec_ec50'):
    cur.execute("ALTER TABLE screen.ba_rec_ec50 Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE screen.ba_rec_ec50 Modify column compound_batch varchar(16)")

    cur.execute("""ALTER TABLE screen.ba_rec_ec50
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
##
if copyTable(engineSCREEN, 'screen.hts_testsets', 'hts_testsets'):
    cur.execute("ALTER TABLE screen.hts_testsets Modify column notebook_ref varchar(16)")

################################################
# GLASS tables

if copyTable(engineGLASS, 'glass.globals', 'globals'):
    try:
        cur.execute("""ALTER TABLE glass.globals CHANGE pkey pkey bigint PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on glass.globals.pkey')
##
if copyTable(engineGLASS, 'glass.hive_stats', 'hive_stats'):
    cur.execute("ALTER TABLE glass.hive_stats Modify column user_id varchar(20)")
    try:
        cur.execute("""ALTER TABLE glass.hive_stats CHANGE pkey pkey bigint PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on glass.hive_stats.pkey')
##
if copyTable(engineGLASS, 'glass.vial', 'vial'):
    cur.execute("ALTER TABLE glass.vial Modify column vial_id varchar(7)")
    cur.execute("ALTER TABLE glass.vial Modify column notebook_ref varchar(30)")
    cur.execute("ALTER TABLE glass.vial Modify column location varchar(40)")
    cur.execute("ALTER TABLE glass.vial Modify column pos varchar(8)")
    try:
        cur.execute("""ALTER TABLE glass.vial CHANGE vial_id vial_id varchar(7) PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on glass.vial.vial_id')
##
if copyTable(engineGLASS, 'glass.vial_log', 'vial_log'):
    cur.execute("ALTER TABLE glass.vial_log Modify column vial_id varchar(7)")
    cur.execute("""ALTER TABLE glass.vial_log
                   ADD FOREIGN KEY (vial_id) REFERENCES glass.vial(vial_id)""")

#
################################################
# SFL_ASSAY tables

if copyTable(engineSFL_ASSAY, 'sfl_assay.rnai', 'rnai'):
    cur.execute("ALTER TABLE sfl_assay.rnai Modify column plate varchar(10)")
    cur.execute("ALTER TABLE sfl_assay.rnai Modify column vendor_reagent_id varchar(20)")

    try:
        cur.execute("""ALTER TABLE sfl_assay.rnai CHANGE pkey pkey bigint PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on sfl_assay.rnai.pkey')

################################################
# LOCTREE tables

if copyTable(engineLOCTREE, 'loctree.location_access', 'location_access'):
    pass
##
if copyTable(engineLOCTREE, 'loctree.location_hierarchy', 'location_hierarchy'):
    try:
        cur.execute("""ALTER TABLE loctree.location_hierarchy CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on loctree.location_hierarchy.pkey')
##
if copyTable(engineLOCTREE, 'loctree.location_type', 'location_type'):
    try:
        cur.execute("""ALTER TABLE loctree.location_type CHANGE type_id type_id int AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on loctree.location_type.type_id')
##
if copyTable(engineLOCTREE, 'loctree.locations', 'locations'):
    cur.execute("ALTER TABLE loctree.locations Modify column parent varchar(40)")
    try:
        cur.execute("""ALTER TABLE loctree.locations CHANGE loc_id loc_id varchar(40) PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on loctree.locations.loc_id')
##
if copyTable(engineLOCTREE, 'loctree.subpositions', 'subpositions'):
    try:
        cur.execute("""CREATE UNIQUE INDEX subpositions_uk_type_seq on loctree.subpositions(type_id, seq)""")
    except Exception as e:
        print(str(e))
        print('Error creating index on loctree.subpositions')

################################################
# MICROTUBE tables

if copyTable(engineMICROTUBE, 'microtube.matrix', 'matrix'):
    cur.execute("ALTER TABLE microtube.matrix Modify column matrix_id varchar(40)")
    cur.execute("ALTER TABLE microtube.matrix Modify column location varchar(80)")

    try:
        cur.execute("""ALTER TABLE microtube.matrix CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on microtube.matrix.pkey')
##
if copyTable(engineMICROTUBE, 'microtube.matrix_position', 'matrix_position'):
    cur.execute("ALTER TABLE microtube.matrix_position Modify column position varchar(10)")

    try:
        cur.execute("""ALTER TABLE microtube.matrix_position CHANGE seq seq int PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on microtube.matrix_position.seq')
##
if copyTable(engineMICROTUBE, 'microtube.matrix_tube', 'matrix_tube'):
    cur.execute("ALTER TABLE microtube.matrix_tube Modify column matrix_id varchar(40)")
    cur.execute("ALTER TABLE microtube.matrix_tube Modify column position varchar(10)")
    cur.execute("ALTER TABLE microtube.matrix_tube Modify column tube_id varchar(10)")

    try:
        cur.execute("""ALTER TABLE microtube.matrix_tube CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on microtube.matrix_tube.pkey')
##
if copyTable(engineMICROTUBE, 'microtube.tube', 'tube'):
    cur.execute("ALTER TABLE microtube.tube Modify column location varchar(60)")
    cur.execute("ALTER TABLE microtube.tube Modify column notebook_ref varchar(30)")
    cur.execute("ALTER TABLE microtube.tube Modify column tube_id varchar(10)")

    try:
        cur.execute("""ALTER TABLE microtube.tube CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on microtube.tube.pkey')
##
if copyTable(engineMICROTUBE, 'microtube.tube_changes', 'tube_changes'):
    cur.execute("ALTER TABLE microtube.tube_changes Modify column location varchar(60)")
    cur.execute("ALTER TABLE microtube.tube_changes Modify column notebook_ref varchar(30)")
    cur.execute("ALTER TABLE microtube.tube_changes Modify column tube_id varchar(10)")

    try:
        cur.execute("""ALTER TABLE microtube.tube_changes CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on microtube.tube_changes.pkey')

################################################
# COOL tables

if copyTable(engineCOOL, 'cool.config', 'config'):
    cur.execute("ALTER TABLE cool.config Modify column config_id varchar(7)")
    cur.execute("ALTER TABLE cool.config Modify column well varchar(10)")
    cur.execute("ALTER TABLE cool.config Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE cool.config Modify column notebook_ref varchar(42)")
##

if copyTable(engineCOOL, 'cool.config_key', 'config_key'):
    cur.execute("ALTER TABLE cool.config_key Modify column config_id varchar(7)")
    #cur.execute("ALTER TABLE cool.config_key Modify column key varchar(32)")
    try:
        cur.execute("""ALTER TABLE cool.config_key CHANGE config_id config_id varchar(7) PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on cool.config_key.config_id')
##

if copyTable(engineCOOL, 'cool.map96to384', 'map96to384'):
    cur.execute("ALTER TABLE cool.map96to384 Modify column well96 varchar(3)")

    cur.execute("ALTER TABLE cool.map96to384 Modify column well384 varchar(4)")
##
if copyTable(engineCOOL, 'cool.plate', 'plate'):
    cur.execute("ALTER TABLE cool.plate Modify column plate_id varchar(7)")
    cur.execute("ALTER TABLE cool.plate Modify column config_id varchar(7)")
    cur.execute("ALTER TABLE cool.plate Modify column location varchar(100)")
    try:
        cur.execute("""ALTER TABLE cool.plate CHANGE plate_id plate_id varchar(7) PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on cool.plate.config_id')
##
if copyTable(engineCOOL, 'cool.plate_type', 'plate_type'):
    try:
        cur.execute("""ALTER TABLE cool.plate_type CHANGE type_id type_id varchar(7) PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on cool.plate_type.type_id')
##
if copyTable(engineCOOL, 'cool.plating_sequence', 'plating_sequence'):
    cur.execute("ALTER TABLE cool.plating_sequence Modify column well varchar(10)")
##
if copyTable(engineCOOL, 'cool.solubility_problem', 'solubility_problem'):
    cur.execute("ALTER TABLE cool.solubility_problem Modify column notebook_ref varchar(20)")
    try:
        cur.execute("""ALTER TABLE cool.solubility_problem CHANGE pkey pkey bigint PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on cool.solubility_problem.pkey')


################################################
# Assay tables

##
if copyTable(engineASSAY, 'assay.lcb_sp', 'lcb_sp'):
    cur.execute("ALTER TABLE assay.lcb_sp Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.lcb_sp Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.lcb_sp Modify column project varchar(40)")
    cur.execute("ALTER TABLE assay.lcb_sp Modify column target varchar(40)")
    cur.execute("ALTER TABLE assay.lcb_sp Modify column operator varchar(40)")
    cur.execute("ALTER TABLE assay.lcb_sp Modify column assay_type varchar(20)")
    cur.execute("ALTER TABLE assay.lcb_sp Modify column eln_id varchar(40)")
    cur.execute("ALTER TABLE assay.lcb_sp Modify column detection_type varchar(20)")

    cur.execute("""create table assay.lcb_sp_invalid as 
    (select * from assay.lcb_sp where compound_batch not in (
    select notebook_ref from bcpvs.batch
    ))""")
    
    cur.execute("""delete from assay.lcb_sp where compound_batch in 
    (select * from
    (select distinct(compound_batch) from assay.lcb_sp where compound_batch not in 
    (select notebook_ref from bcpvs.batch))tmpTbl)
    """)

    cur.execute("""ALTER TABLE assay.lcb_sp
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
##
if copyTable(engineASSAY, 'assay.assay_types', 'assay_types'):
    cur.execute("ALTER TABLE assay.assay_types Modify column assay_type varchar(100)")

    try:
        cur.execute("""ALTER TABLE assay.assay_types CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.assay_types.pkey')
##
if copyTable(engineASSAY, 'assay.ba_bf_camp_dr', 'ba_bf_camp_dr'):
    cur.execute("ALTER TABLE assay.ba_bf_camp_dr Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_dr Modify column compound_id varchar(26)")

    cur.execute("ALTER TABLE assay.ba_bf_camp_dr Modify column project varchar(40)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_dr Modify column method varchar(20)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_dr Modify column target varchar(20)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_dr Modify column operator varchar(20)")
    
    cur.execute("""ALTER TABLE assay.ba_bf_camp_dr
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.ba_bf_camp_dr
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.ba_bf_camp_dr CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.ba_bf_camp_dr.pkey')

##
if copyTable(engineASSAY, 'assay.ba_bf_camp_sp', 'ba_bf_camp_sp'):
    cur.execute("ALTER TABLE assay.ba_bf_camp_sp Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_sp Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_sp Modify column project varchar(20)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_sp Modify column method varchar(20)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_sp Modify column target varchar(20)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_sp Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.ba_bf_camp_sp Modify column notebook_ref varchar(20)")
    
    cur.execute("""ALTER TABLE assay.ba_bf_camp_sp
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.ba_bf_camp_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.ba_bf_camp_sp CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.ba_bf_camp_sp.pkey')
##
if copyTable(engineASSAY, 'assay.ba_biofocus_dr', 'ba_biofocus_dr'):
    cur.execute("ALTER TABLE assay.ba_biofocus_dr Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.ba_biofocus_dr Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.ba_biofocus_dr Modify column project varchar(20)")
    cur.execute("ALTER TABLE assay.ba_biofocus_dr Modify column method varchar(20)")
    cur.execute("ALTER TABLE assay.ba_biofocus_dr Modify column target varchar(20)")
    cur.execute("ALTER TABLE assay.ba_biofocus_dr Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.ba_biofocus_dr Modify column notebook_ref varchar(20)")
    
    cur.execute("""ALTER TABLE assay.ba_biofocus_dr
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.ba_biofocus_dr
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.ba_biofocus_dr CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.ba_biofocus_dr')
##
if copyTable(engineASSAY, 'assay.ba_biofocus_sp', 'ba_biofocus_sp'):
    cur.execute("ALTER TABLE assay.ba_biofocus_sp Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.ba_biofocus_sp Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.ba_biofocus_sp Modify column project varchar(20)")
    cur.execute("ALTER TABLE assay.ba_biofocus_sp Modify column method varchar(20)")
    cur.execute("ALTER TABLE assay.ba_biofocus_sp Modify column target varchar(20)")
    cur.execute("ALTER TABLE assay.ba_biofocus_sp Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.ba_biofocus_sp Modify column notebook_ref varchar(20)")

    cur.execute("""ALTER TABLE assay.ba_biofocus_sp
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.ba_biofocus_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.ba_biofocus_sp CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.ba_biofocus_sp')

##
if copyTable(engineASSAY, 'assay.ba_mds_dr', 'ba_mds_dr'):
    cur.execute("ALTER TABLE assay.ba_mds_dr Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.ba_mds_dr Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.ba_mds_dr Modify column project varchar(20)")
    cur.execute("ALTER TABLE assay.ba_mds_dr Modify column target varchar(80)")
    cur.execute("ALTER TABLE assay.ba_mds_dr Modify column species varchar(20)")
    cur.execute("ALTER TABLE assay.ba_mds_dr Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.ba_mds_dr Modify column notebook_ref varchar(20)")


    cur.execute("""ALTER TABLE assay.ba_mds_dr
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")

    try:
        cur.execute("""ALTER TABLE assay.ba_mds_dr CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.ba_mds_dr')

##
if copyTable(engineASSAY, 'assay.cerep_assays', 'cerep_assays'):
    try:
        cur.execute("""ALTER TABLE assay.cerep_assays CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.cerep_assays')

##
if copyTable(engineASSAY, 'assay.cerep_functional_dr', 'cerep_functional_dr'):
    cur.execute("ALTER TABLE assay.cerep_functional_dr Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.cerep_functional_dr Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.cerep_functional_dr Modify column project varchar(20)")
    cur.execute("ALTER TABLE assay.cerep_functional_dr Modify column target varchar(80)")
    cur.execute("ALTER TABLE assay.cerep_functional_dr Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.cerep_functional_dr Modify column notebook_ref varchar(20)")
    
    cur.execute("""ALTER TABLE assay.cerep_functional_dr
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.cerep_functional_dr
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.cerep_functional_dr CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.cerep_functional_dr')
##
if copyTable(engineASSAY, 'assay.cerep_functional_sp', 'cerep_functional_sp'):
    cur.execute("ALTER TABLE assay.cerep_functional_sp Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.cerep_functional_sp Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.cerep_functional_sp Modify column project varchar(20)")
    cur.execute("ALTER TABLE assay.cerep_functional_sp Modify column target varchar(80)")
    cur.execute("ALTER TABLE assay.cerep_functional_sp Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.cerep_functional_sp Modify column notebook_ref varchar(20)")
    
    cur.execute("""ALTER TABLE assay.cerep_functional_sp
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.cerep_functional_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.cerep_functional_sp CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.cerep_functional_sp')
##
if copyTable(engineASSAY, 'assay.cerep_ki', 'cerep_ki'):
    cur.execute("ALTER TABLE assay.cerep_ki Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.cerep_ki Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.cerep_ki Modify column project varchar(20)")
    cur.execute("ALTER TABLE assay.cerep_ki Modify column target varchar(80)")
    cur.execute("ALTER TABLE assay.cerep_ki Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.cerep_ki Modify column notebook_ref varchar(20)")
    
    cur.execute("""ALTER TABLE assay.cerep_ki
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.cerep_ki
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.cerep_ki CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.cerep_ki')
##
if copyTable(engineASSAY, 'assay.cerep_receptors', 'cerep_receptors'):

    try:
        cur.execute("""ALTER TABLE assay.cerep_receptors CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.cerep_receptors')
##
if copyTable(engineASSAY, 'assay.cerep_screen', 'cerep_screen'):
    cur.execute("ALTER TABLE assay.cerep_screen Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.cerep_screen Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.cerep_screen Modify column project varchar(20)")
    cur.execute("ALTER TABLE assay.cerep_screen Modify column target varchar(80)")
    cur.execute("ALTER TABLE assay.cerep_screen Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.cerep_screen Modify column notebook_ref varchar(20)")
    
    cur.execute("""ALTER TABLE assay.cerep_screen
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.cerep_screen
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.cerep_screen CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.cerep_screen')
##
if copyTable(engineASSAY, 'assay.detection_types', 'detection_types'):

    try:
        cur.execute("""ALTER TABLE assay.detection_types CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.detection_types')
##
if copyTable(engineASSAY, 'assay.lcb_dr', 'lcb_dr'):
    cur.execute("ALTER TABLE assay.lcb_dr Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.lcb_dr Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.lcb_dr Modify column project varchar(40)")
    cur.execute("ALTER TABLE assay.lcb_dr Modify column assay_type varchar(20)")
    cur.execute("ALTER TABLE assay.lcb_dr Modify column detection_type varchar(20)")
    cur.execute("ALTER TABLE assay.lcb_dr Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.lcb_dr Modify column eln_id varchar(20)")
    cur.execute("""create table assay.lcb_dr_invalid as 
    (select * from assay.lcb_dr where compound_batch not in (
    select notebook_ref from bcpvs.batch
    ))""")
    
    cur.execute("""delete from assay.lcb_dr where compound_batch in 
    (select * from
    (select distinct(compound_batch) from assay.lcb_dr where compound_batch not in 
    (select notebook_ref from bcpvs.batch))tmpTbl)
    """)
    
    cur.execute("""ALTER TABLE assay.lcb_dr
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.lcb_dr
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.lcb_dr CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.lcb_dr')
##
if copyTable(engineASSAY, 'assay.lcb_sp_schneider_cysm', 'lcb_sp_schneider_cysm'):
    cur.execute("ALTER TABLE assay.lcb_sp_schneider_cysm Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.lcb_sp_schneider_cysm Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.lcb_sp_schneider_cysm Modify column project varchar(40)")
    cur.execute("ALTER TABLE assay.lcb_sp_schneider_cysm Modify column target varchar(80)")
    cur.execute("ALTER TABLE assay.lcb_sp_schneider_cysm Modify column assay_type varchar(20)")
    cur.execute("ALTER TABLE assay.lcb_sp_schneider_cysm Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.lcb_sp_schneider_cysm Modify column detection_type varchar(20)")
    cur.execute("ALTER TABLE assay.lcb_sp_schneider_cysm Modify column eln_id varchar(20)")
    
    cur.execute("""create table assay.lcb_sp_schneider_cysm_invalid as 
    (select * from assay.lcb_sp_schneider_cysm where compound_batch not in (
    select notebook_ref from bcpvs.batch
    ))""")
    
    cur.execute("""delete from assay.lcb_sp_schneider_cysm where compound_batch in 
    (select * from
    (select distinct(compound_batch) from assay.lcb_sp_schneider_cysm where compound_batch not in 
    (select notebook_ref from bcpvs.batch))tmpTbl)
    """)

    cur.execute("""ALTER TABLE assay.lcb_sp_schneider_cysm
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")

    try:
        cur.execute("""ALTER TABLE assay.lcb_sp_schneider_cysm CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.lcb_sp_schneider_cysm')

##
if copyTable(engineASSAY, 'assay.metabolite', 'metabolite'):
    cur.execute("ALTER TABLE assay.metabolite Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.metabolite Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.metabolite Modify column project varchar(40)")
    cur.execute("ALTER TABLE assay.metabolite Modify column operator varchar(20)")
    cur.execute("ALTER TABLE assay.metabolite Modify column notebook_ref varchar(20)")
    
    cur.execute("""ALTER TABLE assay.metabolite
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.metabolite
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.metabolite CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.metabolite')
##
if copyTable(engineASSAY, 'assay.upstate_sp', 'upstate_sp'):
    cur.execute("ALTER TABLE assay.upstate_sp Modify column compound_batch varchar(16)")
    cur.execute("ALTER TABLE assay.upstate_sp Modify column compound_id varchar(26)")
    cur.execute("ALTER TABLE assay.upstate_sp Modify column project varchar(40)")
    cur.execute("ALTER TABLE assay.upstate_sp Modify column target varchar(80)")
    cur.execute("ALTER TABLE assay.upstate_sp Modify column method varchar(20)")
    cur.execute("ALTER TABLE assay.upstate_sp Modify column operator varchar(20)")

    cur.execute("""ALTER TABLE assay.upstate_sp
                   ADD FOREIGN KEY (compound_batch) REFERENCES bcpvs.batch(notebook_ref)""")
    cur.execute("""ALTER TABLE assay.upstate_sp
                   ADD FOREIGN KEY (compound_id) REFERENCES bcpvs.compound(compound_id)""")

    try:
        cur.execute("""ALTER TABLE assay.upstate_sp CHANGE pkey pkey bigint AUTO_INCREMENT PRIMARY KEY""")
    except Exception as e:
        print(str(e))
        print('Error creating index on assay.upstate_sp')

################################################
# Chemspec tables

if copyTable(engineCHEMSPEC, 'chemspec.chem_info', 'chem_info'):
    cur.execute("ALTER TABLE chemspec.chem_info Modify column regno bigint")
    cur.execute("ALTER TABLE chemspec.chem_info Modify column jpage varchar(20)")
    cur.execute("ALTER TABLE chemspec.chem_info Modify column solvent varchar(50)")
    cur.execute("ALTER TABLE chemspec.chem_info Modify column project varchar(40)")
    cur.execute("ALTER TABLE chemspec.chem_info Modify column library_id varchar(80)")
    
    try:
        cur.execute("""CREATE UNIQUE INDEX regno_pk ON chemspec.chem_info(regno)""")
    except Exception as e:
        print(str(e))
        print('Error creating index on chemspec.chem_info.regno')

##
if copyTable(engineCHEMSPEC, 'chemspec.analytical_sum', 'analytical_sum'):
    cur.execute("ALTER TABLE chemspec.analytical_sum Modify column regno bigint")

    cur.execute("""delete from chemspec.analytical_sum where regno in 
    (select * from
    (select distinct(regno) from chemspec.analytical_sum where regno not in 
    (select regno from chemspec.chem_info))tmpTbl)
    """)

    cur.execute("""ALTER TABLE chemspec.analytical_sum
                   ADD FOREIGN KEY (regno) REFERENCES chemspec.chem_info(regno)""")

##
if copyTable(engineCHEMSPEC, 'chemspec.chemspec_salt', 'chemspec_salt'):
    cur.execute("ALTER TABLE chemspec.chemspec_salt Modify column suffix varchar(10)")

##
if copyTable(engineCHEMSPEC, 'chemspec.chns_log', 'chns_log'):
    cur.execute("ALTER TABLE chemspec.chns_log Modify column regno bigint")
    cur.execute("""ALTER TABLE chemspec.chns_log
                   ADD FOREIGN KEY (regno) REFERENCES chemspec.chem_info(regno)""")

##
if copyTable(engineCHEMSPEC, 'chemspec.color_tbl', 'color_tbl'):
    cur.execute("ALTER TABLE chemspec.color_tbl Modify column color varchar(20)")

##
if copyTable(engineCHEMSPEC, 'chemspec.experiment', 'experiment'):
    cur.execute("ALTER TABLE chemspec.experiment Modify column technique varchar(20)")

##
if copyTable(engineCHEMSPEC, 'chemspec.ftir_log', 'ftir_log'):
    cur.execute("ALTER TABLE chemspec.ftir_log Modify column regno bigint")
    cur.execute("""ALTER TABLE chemspec.ftir_log
                   ADD FOREIGN KEY (regno) REFERENCES chemspec.chem_info(regno)""")

##
if copyTable(engineCHEMSPEC, 'chemspec.hplc_uv_log', 'hplc_uv_log'):
    cur.execute("ALTER TABLE chemspec.hplc_uv_log Modify column regno bigint")

    cur.execute("""delete from chemspec.hplc_uv_log where regno in 
    (select * from
    (select distinct(regno) from chemspec.hplc_uv_log where regno not in 
    (select regno from chemspec.chem_info))tmpTbl)
    """)
    cur.execute("""ALTER TABLE chemspec.hplc_uv_log
                   ADD FOREIGN KEY (regno) REFERENCES chemspec.chem_info(regno)""")

##
if copyTable(engineCHEMSPEC, 'chemspec.ms_log', 'ms_log'):
    cur.execute("ALTER TABLE chemspec.ms_log Modify column regno bigint")
    cur.execute("""ALTER TABLE chemspec.ms_log
                   ADD FOREIGN KEY (regno) REFERENCES chemspec.chem_info(regno)""")

##
if copyTable(engineCHEMSPEC, 'chemspec.nmr_log', 'nmr_log'):
    cur.execute("ALTER TABLE chemspec.nmr_log Modify column regno bigint")
    cur.execute("""delete from chemspec.nmr_log where regno in 
    (select * from
    (select distinct(regno) from chemspec.nmr_log where regno not in 
    (select regno from chemspec.chem_info))tmpTbl)
    """)
    cur.execute("""ALTER TABLE chemspec.nmr_log
                   ADD FOREIGN KEY (regno) REFERENCES chemspec.chem_info(regno)""")

##
if copyTable(engineCHEMSPEC, 'chemspec.other_log', 'other_log'):
    cur.execute("ALTER TABLE chemspec.other_log Modify column regno bigint")
    cur.execute("""ALTER TABLE chemspec.other_log
                   ADD FOREIGN KEY (regno) REFERENCES chemspec.chem_info(regno)""")

##
if copyTable(engineCHEMSPEC, 'chemspec.solvent_tbl', 'solvent_tbl'):
    cur.execute("ALTER TABLE chemspec.solvent_tbl Modify column solvent varchar(20)")

##
if copyTable(engineCHEMSPEC, 'chemspec.spec_analysis', 'spec_analysis'):
    cur.execute("ALTER TABLE chemspec.spec_analysis Modify column regno bigint")
    cur.execute("""ALTER TABLE chemspec.spec_analysis
                   ADD FOREIGN KEY (regno) REFERENCES chemspec.chem_info(regno)""")

##
if copyTable(engineCHEMSPEC, 'chemspec.storage_tbl', 'storage_tbl'):
    cur.execute("ALTER TABLE chemspec.storage_tbl Modify column storage varchar(30)")

################################################
# Close
pt.close_connection(con)

quit()


