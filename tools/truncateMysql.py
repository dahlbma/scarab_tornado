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


if allTabs:
    truncateTable("bcpvs.batch")
    truncateTable("bcpvs.batch_invalid")
    truncateTable("bcpvs.compound")


cur.execute('SET FOREIGN_KEY_CHECKS = 1')
