import MySQLdb 
import config
import logging

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')


def setup_logger(name, log_file, level=logging.INFO):
    handler = logging.FileHandler(log_file)        
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

scarabLogger = setup_logger('scara_logger', 'chemRegSqlError.txt')

class DisconnectSafeCursor(object):
    db = None
    cursor = None

    def __init__(self, db, cursor):
        self.db = db
        self.cursor = cursor

    def close(self):
        self.cursor.close()


    def ping(self, *args, **kwargs):
        ret = ''
        try:
            # Pass False (or nothing) to avoid the deprecated auto-reconnect
            self.db.conn.ping(False)
        except (MySQLdb.OperationalError, MySQLdb.InterfaceError) as e:
            scarabLogger.error(f"Connection lost, attempting manual reconnect: {str(e)}")
            try:
                self.db.reconnect()
                # Update the cursors after reconnecting
                self.cursor = self.db.cur
                scarabLogger.info("Reconnection successful.")
            except Exception as re_err:
                scarabLogger.error(f"Reconnection failed: {str(re_err)}")
                ret = 'error'
        return ret

    
    def execute(self, *args, **kwargs):
        try:
            sSql = args[0].encode('utf-8', 'replace').decode('utf-8')
        except:
            sSql = args[0]
        try:
            return self.cursor.execute(*args, **kwargs)
        except Exception as e:
            scarabLogger.error(str(e))
            scarabLogger.error(args)
            return -1
        
    def fetchone(self):
        return self.cursor.fetchone()

    def fetchall(self):
        return self.cursor.fetchall()

    def description(self):
        return self.cursor.description


class DisconnectSafeConnection(object):
    connect_args = None
    connect_kwargs = None
    conn = None
    cur = None
    
    def __init__(self, *args, **kwargs):
        self.connect_args = args
        self.connect_kwargs = kwargs
        self.reconnect()

    def reconnect(self):
        try:
            self.conn.close()
        except:
            pass
        
        self.conn = MySQLdb.connect(
            host=config.database['host'],
            user=config.database['user'],
            passwd=config.database['password'],
            database=config.database['db'],
            charset='utf8mb4',  # Use utf8mb4 charset
            use_unicode=True    # Use Unicode
        )
        self.conn.autocommit(True)


    def cursor(self, *args, **kwargs):
        self.cur = self.conn.cursor(*args, **kwargs)
        return DisconnectSafeCursor(self, self.cur)

    def commit(self):
        self.conn.commit()

    def rollback(self):
        self.conn.rollback()

disconnectSafeConnect = DisconnectSafeConnection
