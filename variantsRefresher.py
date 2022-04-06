#!/usr/local/bin/python3

import re, json, yaml, sys, datetime, statistics
from os import path, environ, pardir
from pymongo import MongoClient
from isodate import date_isoformat
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
    set_processing_modes(byc)

    select_dataset_ids(byc)
    check_dataset_ids(byc)

    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d ...")
        exit()

    if byc["test_mode"]:
        max_count = 10
    else:
        max_count = int(byc["args"].randno)

    mongo_client = MongoClient( )

    for ds_id in byc["dataset_ids"]:

        data_db = mongo_client[ ds_id ]
        var_coll = data_db[ "variants" ]
        bios_coll = data_db[ "biosamples" ]
        no =  bios_coll.estimated_document_count()
        if max_count < 1:
            max_count = no
        count = 0
        bar = Bar("{} vars".format(ds_id), max = no, suffix='%(percent)d%%'+" of {} biosamples".format(no) )
        for bios in bios_coll.find({}):
            bios_id = bios["id"]
            ind_id = bios.get("individual_id", "")

            for v_p in var_coll.find({"biosample_id": bios_id}):
                var_coll.update_one( { "_id": v_p["_id"] }, { '$set': {"individual_id": ind_id} }  )            
            bar.next()
        bar.finish()

    exit()

    v_d = byc["variant_definitions"]

    for ds_id in byc["dataset_ids"]:

        v_short = 0
        v_no_type = 0

        data_db = mongo_client[ ds_id ]
        var_coll = data_db[ "variants" ]
        bios_coll = data_db[ "biosamples" ]
        no =  var_coll.estimated_document_count()
        if max_count < 1:
            max_count = no
        count = 0
        bar = Bar("{} vars".format(ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )
        for v_p in var_coll.find({}):

            if count > max_count:
                break
            count += 1

            update_obj = vrsify_variant(v_p, v_d)
            update_obj.update({ "updated": datetime.datetime.now().isoformat() })
            
            # for p in ["callset_id", "biosample_id", "variant_internal_id", "variant_type", "variant_state", "reference_bases", "alternate_bases", "info"]:
            #     p_v = v_p.get(p, None)
            #     if p_v is not None:
            #         update_obj.update({p:p_v})

            if not "start" in v_p:
                break

            if not byc["test_mode"]:
                var_coll.update_one( { "_id": v_p["_id"] }, { '$set': update_obj }  )
            else:
                prjsonnice(update_obj)
            
            bar.next()
        bar.finish()

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
