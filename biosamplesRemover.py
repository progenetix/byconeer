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

## `biosamplesRemover`

* `biosamplesRefresher.py -d progenetix -m stats`

"""

################################################################################
################################################################################
################################################################################

def main():

    biosamples_remover()

################################################################################

def biosamples_remover():

    initialize_service(byc, "biosamples_refresher")
    get_args(byc)

    if not byc["args"].query:
        print('¡¡¡ No sample query specified - use `-q {"biosamples":{ ...} }"')
        exit()

    query = json.loads(byc["args"].query)

    set_processing_modes(byc)

    select_dataset_ids(byc)
    check_dataset_ids(byc)
    parse_filters(byc)
    parse_variants(byc)
    initialize_beacon_queries(byc)
    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d ...")
        exit()
    if len(byc["dataset_ids"]) > 1:
        print("Only one dataset should be provided with -d ...")
        exit()

    ds_id =  byc["dataset_ids"][0]
    junk_id = byc["config"]["dataset_junk"]

    counts = {
        "biosamples": 0,
        "callsets": 0,
        "variants": 0,
        "individuals": 0
    }

    bios_query = {}
    if "biosamples" in query:
        bios_query = query["biosamples"]

    print(bios_query)

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    bios_coll = data_db[ "biosamples" ]
    cs_coll = data_db[ "callsets" ]
    v_coll = data_db[ "variants" ]
    ind_coll = data_db[ "individuals" ]

    bs_uids = []

    for bs in bios_coll.find (bios_query, {"_id":1, "id":1} ):
        bs_uids.append(bs["_id"])

    bios_no =  len(bs_uids)

    if bios_no < 1:
        print('!!! No biosamples have been found ¯\_(ツ)_/¯ ')
        exit()

    print("??? {} biosamples have been found => move to **{}**?".format(bios_no, junk_id))


    reason = input('Please type a brief note for the reason:'+"\n=> ")
    if len(reason) < 5:
        print("Not enough ... Try again!")
        exit()

    resp_test = "JUNK {} SAMPLES".format(bios_no)
    resp = input("Please type \"{}\" or otherwise something else to abort:\n=> ".format(resp_test))
    if not resp_test in resp:
        exit()

    print("... now collecting data & junking it ...")

    # individuals have to be collected and checked - they may have other biosamples
    ind_ids = []

    junk_db = data_client[ junk_id ]
    bios_jcoll = junk_db[ "biosamples" ]
    cs_jcoll = junk_db[ "callsets" ]
    v_jcoll = junk_db[ "variants" ]
    ind_jcoll = junk_db[ "individuals" ]

    j_info = { "junk_note": reason, "junk_date": date_isoformat(datetime.datetime.now()) }


    for bsuid in bs_uids:
        s = bios_coll.find_one({ "_id":bsuid })
        ind_ids.append(s["individual_id"])
        counts["biosamples"] += 1

        s = _update_info(s, j_info)

        if not byc["test_mode"]:
            bios_jcoll.insert_one(s)
            bios_coll.delete_one({"_id":bsuid})

        for cs in cs_coll.find( {"biosample_id":s["id"] } ):
            cs = _update_info(cs, j_info)
            counts["callsets"] += 1
            if not byc["test_mode"]:
                cs_jcoll.insert_one(cs)
                cs_coll.delete_one({"_id":cs["_id"]})

        for v in v_coll.find({"biosample_id":s["id"]}):
            v = _update_info(v, j_info)
            counts["variants"] += 1
            if not byc["test_mode"]:
                v_jcoll.insert_one(v)
                v_coll.delete_one({"_id":v["_id"]})

    for ind_id in ind_ids:

        bios = bios_coll.find_one({ "individual_id":ind_id })

        if not bios:

            ind = ind_coll.find_one({"id":ind_id})
            ind = _update_info(ind, j_info)
            counts["individuals"] +=1
            if not byc["test_mode"]:
                ind_jcoll.insert_one(ind)
                ind_coll.delete_one({"_id":ind["_id"]})


    for k, n in counts.items():
        print("=> junked {}: {}".format(k, n))

    print("=> The documents have been moved to th ecorresponding collections in {}".format(junk_id))

################################################################################

def _update_info(o, j_info):

    i = o.get("info", {})
    for k, v in j_info.items():
        i.update({k:v})
    o.update({"info": i})
    return o


################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
