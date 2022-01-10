#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from isodate import date_isoformat
from pymongo import MongoClient
import argparse
import statistics
from progress.bar import Bar

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""

## `biosamplesRefresher`

"""

################################################################################
################################################################################
################################################################################

def main():

    variants_refresher()

################################################################################

def variants_refresher():

    initialize_service(byc)
    get_args(byc)
    set_test_mode(byc)

    select_dataset_ids(byc)
    check_dataset_ids(byc)

    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d ...")
        exit()

    mongo_client = MongoClient( )

    min_l = byc["this_config"]["refreshing"]["cnv_min_length"]

    for ds_id in byc["dataset_ids"]:

        v_short = 0
        v_no_type = 0

        data_db = mongo_client[ ds_id ]
        var_coll = data_db[ "variants" ]
        no =  var_coll.estimated_document_count()
        bar = Bar("{} vars".format(ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )
        for v in var_coll.find({}):
            update_obj = { "id": str(v["_id"]) }

            if not byc["test_mode"]:
                var_coll.update_one( { "_id": v["_id"] }, { '$set': update_obj }  )
            
            bar.next()
        bar.finish()

        print("{} {} variants had no type {}".format(v_no_type, ds_id, min_l))
        print("{} {} CNV variants had a length below {}".format(v_short, ds_id, min_l))

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
