#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from isodate import date_isoformat
from pymongo import MongoClient
import argparse
from progress.bar import Bar
import time

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""

## `frequencymapsCreator`

"""

################################################################################
################################################################################
################################################################################

def main():
    frequencymaps_creator()

################################################################################

def frequencymaps_creator():

    initialize_service(byc)
    run_beacon_init_stack(byc)

    if len(byc["dataset_ids"]) > 1:
        print("Please give only one dataset using -d")
        exit()

    ds_id = byc["dataset_ids"][0]

    # for non-standard CNV binning
    genome_binning_from_args(byc)
    generate_genomic_intervals(byc)
    
    print("=> Using data values from {}".format(ds_id))

    coll_client = MongoClient()
    coll_coll = coll_client[ ds_id ][ byc["config"]["collations_coll"] ]

    fm_client = MongoClient()
    fm_coll = fm_client[ ds_id ][ byc["config"]["frequencymaps_coll"] ]

    bios_client = MongoClient()
    bios_coll = bios_client[ ds_id ][ byc["config"]["collations_source"] ]

    ind_client = MongoClient()
    ind_coll = ind_client[ ds_id ]["individuals"]

    data_client = MongoClient()
    cs_coll = data_client[ ds_id ]["callsets"]

    id_query = {}
    id_ql = []

    if len(byc["filters"]) > 0:
        f_l = []
        for c_t in byc["filters"]:
            f_l.append( { "id": { "$regex": "^"+c_t["id"] } })
        if len(f_l) > 1:
            id_ql.append( { "$or": f_l } )
        else:
            id_ql.append(f_l[0])

    if len(id_ql) == 1:
        id_query = id_ql[0]
    elif len(id_ql) > 1:
        id_query = { "$and":id_ql }

    coll_ids = coll_coll.distinct("_id", id_query)
    coll_no = len(coll_ids)
   
    print("Writing {} {} fMaps for {} intervals".format(coll_no, ds_id, len(byc["genomic_intervals"])))

    coll_i = 0

    for c_id in coll_ids:

        coll = coll_coll.find_one({"_id":c_id})

        if not coll:
            print("?????? some error - collation {} not found !!!".format(c_id))
            continue

        pre, code = re.split("[:-]", coll["id"], 1)
        coll_type = coll.get("collation_type", "undefined")

        exclude_normals = True

        for normal in byc["service_config"]["keep_normals"]:
            if normal in coll["id"]:
                print("---> keeping normals for {}".format(coll["id"]))
                exclude_normals = False

        db_key = coll["db_key"]

        coll_i += 1

        query = { db_key: { '$in': coll["child_terms"] } }
        bios_no, cs_cursor = _cs_cursor_from_bios_query(bios_coll, ind_coll, cs_coll, coll["id"], coll["scope"], query, exclude_normals)
        cs_no = len(list(cs_cursor))

        if cs_no < 1:
            continue

        i_t = coll_i % 100
        start_time = time.time()
        if i_t == 0 or cs_no > 1000:
            print("{}: {} bios, {} cs\t{}/{}\t{:.1f}%".format(coll["id"], bios_no, cs_no, coll_i, coll_no, 100*coll_i/coll_no))

        update_obj = {
            "id": coll["id"],
            "label": coll["label"],
            "dataset_id": coll["dataset_id"],
            "scope": coll["scope"],
            "db_key": coll["db_key"],
            "collation_type": coll["collation_type"],
            "child_terms": coll["child_terms"],
            "updated": date_isoformat(datetime.datetime.now()),
            "counts": {"biosamples": bios_no, "callsets": cs_no },
            "frequencymap": {
                "interval_count": len(byc["genomic_intervals"]),
                "binning": byc["genome_binning"],
                "biosample_count": bios_no,
                "analysis_count": cs_no,
                "intervals": interval_counts_from_callsets(cs_cursor, byc)
            }
        }

        proc_time = time.time() - start_time
        if cs_no > 1000:
            print(" => Processed in {:.2f}s: {:.4f}s per callset".format(proc_time, (proc_time/cs_no)))

        if not byc["test_mode"]:
            print("Updating {}...".format(coll["id"]))
            fm_coll.delete_one( { "id": coll["id"] } )
            fm_coll.insert_one( update_obj )

        if coll["code_matches"] > 0:
            if cs_no != coll["code_matches"]:
                query_cm = { db_key: coll["id"] }
                bios_no_cm, cs_cursor_cm = _cs_cursor_from_bios_query(bios_coll, ind_coll, cs_coll, coll["id"], coll["scope"], query_cm)
                cs_no_cm = len(list(cs_cursor_cm))
                if cs_no_cm > 0:

                    cm_obj = { "frequencymap_codematches": {
                            "interval_count": len(byc["genomic_intervals"]),
                            "binning": byc["genome_binning"],
                            "biosample_count": bios_no_cm,
                            "analysis_count": cs_no_cm,
                            "intervals": interval_counts_from_callsets(cs_cursor_cm, byc)
                        }
                    }

                    print("{}: {} exact of {} total code matches".format(coll["id"], cs_no_cm, cs_no))

                    if not byc["test_mode"]:
                        fm_coll.update_one( { "id": coll["id"] }, { '$set': cm_obj }, upsert=False )

################################################################################

def _cs_cursor_from_bios_query(bios_coll, ind_coll, cs_coll, coll_id, scope, query, exclude_normals=True):

    if scope == "individuals":
        ind_ids = ind_coll.distinct( "id" , query )
        bios_ids = bios_coll.distinct( "id" , {"individual_id":{"$in": ind_ids } } )
    elif scope == "callsets":
        bios_ids = cs_coll.distinct( "biosample_id" , query )
    else:
        bios_ids = bios_coll.distinct( "id" , query )

    pre_b = len(bios_ids)

    # for most entities samples labeled as "normal" will be excluded for frequency calculations
    if exclude_normals:
        bios_ids = bios_coll.distinct( "id" , { "id": { "$in": bios_ids } , "biosample_status.id": {"$ne": "EFO:0009654" }} )
    bios_no = len(bios_ids)
    
    if pre_b > bios_no:
        print("WARNING: {} samples for {}, while {} after excluding normals by EFO:0009654".format(pre_b, coll_id, bios_no))
       
    cs_query = { "biosample_id": { "$in": bios_ids } }
    cs_cursor = cs_coll.find(cs_query)

    return bios_no, cs_cursor

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
