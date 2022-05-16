#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
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

    initialize_service(byc, "callsets")
    run_beacon_init_stack(byc)

    if len(byc["dataset_ids"]) > 1:
        print("Please give only one dataset using -d")
        exit()

    ds_id = byc["dataset_ids"][0]

    # for non-standard CNV binning
    genome_binning_from_args(byc)
    generate_genomic_intervals(byc)
    
    print("=> Using data values from {}".format(ds_id))

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    cs_coll = data_db[ "callsets" ]
    v_coll = data_db[ "variants" ]

    execute_bycon_queries( ds_id, byc )

    ds_results = byc["dataset_results"][ds_id]

    if not "callsets._id" in ds_results.keys():
        cs_ids = cs_coll.distinct("_id", {})
        print("¡¡¡ Using all {} callsets from {} !!!".format(len(cs_ids), ds_id))
    else:
        cs_ids = ds_results["callsets._id"]["target_values"]

    print("Re-generating statusmaps with {} intervals for {} callsets...".format(len(byc["genomic_intervals"]), len(cs_ids)))

    no =  len(cs_ids)
    bar = Bar("{} callsets".format(ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )

    for _id in cs_ids:

        cs = cs_coll.find_one( { "_id": _id } )
        csid = cs["id"]

        bar.next()

        # only the defined parameters will be overwritten
        cs_update_obj = { "info": cs.get("info", {}) }
        cs["info"].pop("statusmaps", None)
        cs["info"].pop("cnvstatistics", None)

        maps, cs_cnv_stats, cs_chro_stats = interval_cnv_arrays(v_coll, { "callset_id": csid }, byc)
        cs_update_obj.update({"cnv_statusmaps": maps})
        cs_update_obj.update({"cnv_stats": cs_cnv_stats})
        cs_update_obj.update({"cnv_chro_stats": cs_chro_stats})
        cs_update_obj.update({ "updated": datetime.datetime.now().isoformat() })

        if not byc["test_mode"]:
            cs_coll.update_one( { "_id": _id }, { '$set': cs_update_obj }  )
        else:
            prjsonnice(maps)

        ####################################################################
        ####################################################################
        ####################################################################

    bar.finish()

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
