import datetime

################################################################################

def set_db_object_defaults(dbobj, scope, byc):

    defaults = byc["db_object_defaults"][scope]
    for k, v in defaults:
        dbobj.update({k:v})
    dbobj.update({ "updated": datetime.datetime.now().isoformat() })

    return dbobj

################################################################################
