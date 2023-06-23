#!/usr/local/bin/python3

import re, json, yaml, sys, datetime, statistics
from os import path, environ, pardir
from pymongo import MongoClient
from isodate import date_isoformat
from progress.bar import Bar

from bycon import *

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

    initialize_bycon_service(byc)
    parse_variant_parameters(byc)
    select_dataset_ids(byc)
    
    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d ...")
        exit()

    if byc["test_mode"]:
        max_count = 10000000
    else:
        max_count = int(byc["args"].randno)

    mongo_client = MongoClient( )

    for ds_id in byc["dataset_ids"]:

        v_short = 0
        v_no_type = 0

        data_db = mongo_client[ ds_id ]
        var_coll = data_db[ "variants" ]

        no =  var_coll.estimated_document_count()
        if max_count < 1:
            max_count = no
        count = 0
        bar = Bar("{} vars".format(ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )
        for v_p in var_coll.find({}):

            if count > max_count:
                break
            count += 1

            update_obj = v23(v_p)
            update_obj.update({
                "updated": datetime.datetime.now().isoformat()
            })

            if not byc["test_mode"]:
                var_coll.update_one( { "_id": v_p["_id"] }, { '$set': update_obj }  )
            # else:
            #     prjsonnice(update_obj)
            
            bar.next()
        bar.finish()

################################################################################

################################################################################

def v23(v):

    v_d = byc["variant_definitions"]

    if "M:" in v.get("variant_internal_id", ""):
        v.update({"variant_internal_id": re.sub("M:", "MT:", v["variant_internal_id"])})

    # protection against re-doing the conversion (which would lose data)
    if "info" in v:
        if "v23" in v["info"].get("version", ""):
            return v
    else:
        v.update({"info":{}})


    if "MT:" in v.get("variant_internal_id", ""):
        v["location"].update({"sequence_id": "refseq:NC_012920.1"})

    if v["location"]["sequence_id"] is None:
        prjsonnice(v)
        exit()

    v.update({
        "location": {
            "sequence_id": v["location"]["sequence_id"],
            "chromosome": v_d["refseq_chronames"][ v["location"]["sequence_id"] ],
            "start": int(v["location"]["interval"]["start"]["value"]),
            "end": int(v["location"]["interval"]["end"]["value"])
        },
    })

    v["info"].update({"version":"v23"})

    if "state" in v:
        s_t = v["state"].get("type", "")
        if "LiteralSequenceExpression" in s_t:
            s_s = v["state"].get("sequence", "")
            v.update({
                "variant_state": {
                    "id": "SO:0001059",
                    "label": "sequence_alteration"
                },
                "sequence": s_s,
                "reference_sequence": v.get("reference_bases")
            })
            v.pop("reference_bases", None)

    return v

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
