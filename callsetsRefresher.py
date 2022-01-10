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

################################################################################
"""
* `byconeer/callsetsRefresher.py -d 1000genomesDRAGEN -s variants`
* `byconeer/callsetsRefresher.py -d progenetix -s biosamples -f "icdom-81703"`
* `byconeer/callsetsRefresher.py`
  - default; new statusmaps for all `progenetix` callsets
"""
################################################################################

def main():

    callsets_refresher()

################################################################################

def callsets_refresher():

    initialize_service(byc)
    get_args(byc)
    set_test_mode(byc)

    select_dataset_ids(byc)
    check_dataset_ids(byc)
    parse_filters(byc)
    parse_variants(byc)
    initialize_beacon_queries(byc)

    generate_genomic_intervals(byc)

    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d => using progenetix")

    generate_genomic_intervals(byc)

    for ds_id in byc["dataset_ids"]:
        print("=> Using callset_id values from {}.{}".format(ds_id, byc["args"].source))
        _process_dataset(ds_id, byc)

################################################################################
################################################################################
################################################################################

def _process_dataset(ds_id, byc):

    no_cs_no = 0
    no_stats_no = 0

    cs_id_ks = {}

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    cs_coll = data_db[ "callsets" ]
    v_coll = data_db[ "variants" ]

    if byc["args"].source == "biosamples":
        bios_query = {}
        if "biosamples" in byc["queries"]:
            bios_query = byc["queries"]["biosamples"]
        bios_coll = data_db[ "biosamples" ]
        bs_ids = []
        for bs in bios_coll.find (bios_query, {"id":1} ):
            bs_ids.append(bs["id"])

        for bsid in bs_ids:
            cs_query = { "biosample_id": bsid }
            for cs in cs_coll.find (cs_query ):
                cs_id_ks.update({cs["callset_id"]: cs["biosample_id"]})
    elif byc["args"].source == "variants":
        v_query = {}
        if "variants" in byc["queries"]:
            v_query = byc["queries"]["variants"]
        for v in v_coll.find (v_query ):
            cs_id_ks.update({v["callset_id"]: v["biosample_id"]})
    else:
        for cs in cs_coll.find({}):
            cs_id_ks.update({cs["id"]: cs["biosample_id"]})

    cs_ids = list(cs_id_ks.keys())

    no =  len(cs_ids)

    bar = Bar("{} callsets from {}".format(no, ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )
    
    for csid in cs_ids:

        cs = cs_coll.find_one( { "id": csid } )

        bar.next()

        if not cs:
            no_cs_no += 1
            cs_update_obj = {
                "id": csid,
                "biosample_id": cs_id_ks[csid],
                "info": {}
            }
        elif not "info" in cs:
            cs_update_obj = { "info":{} }
        else:
            cs_update_obj = { "info": cs["info"] }

        maps, cs_cnv_stats = interval_cnv_arrays(v_coll, { "callset_id": csid }, byc)

        cs_update_obj["info"].update({"statusmaps": maps})
        cs_update_obj["info"].update({"cnvstatistics": cs_cnv_stats})
        cs_update_obj.update({ "updated": date_isoformat(datetime.datetime.now()) })

        if not byc["test_mode"]:
            if not cs:
                cs_coll.insert_one( cs_update_obj  )
            else:
                cs_coll.update_one( { "_id": cs["_id"] }, { '$set': cs_update_obj }  )
        else:
            print(json.dumps(camelize(maps), sort_keys=True, default=str))

        ####################################################################
        ####################################################################
        ####################################################################

    bar.finish()

    print("{} {} biosamples had no callsets".format(no_cs_no, ds_id))

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
