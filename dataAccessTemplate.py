#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from pymongo import MongoClient
import statistics
from progress.bar import Bar

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

def main():

    data_retriever()

################################################################################

# ./dataAccessTemplate.py --filters "pgx:icdom-80702"

def data_retriever():

    """podmd
    One can specify a separate `../config/data_retriever.yaml` config file
    which would be automatically read in:
    * all parameters would be available under `byc["service_config"]`
    * parameters under the `defaults` root parameter would be assigned into
      `byc[ ... ]`; e.g. `byc["response_entity"]` could be defined there
    """

    initialize_service(byc, "biosamples")
    run_beacon_init_stack(byc)

    ds_id = set_single_dataset(byc)

    # for non-standard CNV binning
    genome_binning_from_args(byc)
    generate_genomic_intervals(byc)
    
    print("=> Using data values from {}".format(ds_id))

    execute_bycon_queries( ds_id, byc )
    ds_results = byc["dataset_results"][ds_id]

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    bios_coll = data_db[ "biosamples" ]
    ind_coll = data_db[ "individuals" ]
    cs_coll = data_db[ "callsets" ]
    v_coll = data_db[ "variants" ]

    for k in ds_results.keys():
        if ds_results[k]["target_count"] > 0:
            print("==> available data access key: {}".format(k))

    if not "biosamples._id" in ds_results.keys():
        bios_ids = bios_coll.distinct("_id", {})
        print("¡¡¡ Using all {} callsets from {} !!!".format(len(bios_ids), ds_id))
    else:
        bios_ids = ds_results["biosamples._id"]["target_values"]

    bios_no = len(bios_ids)
    print("Processing data with {} intervals for {} samples...".format(len(byc["genomic_intervals"]), bios_no))

    ############################################################################
    # from her on build your script ...
    ############################################################################

    bar = Bar("{} samples".format(ds_id), max = bios_no, suffix='%(percent)d%%'+" of "+str(bios_no) )

    # if other datya access keys had been initiated through a search one could
    # use the values from there, e.g. ds_results["callsets._id"]["target_values"]
    # for accessing callsets...

    for _id in bios_ids:

        bios = bios_coll.find_one( { "_id": _id } )


        bar.next()

        # only the defined parameters will be overwritten
        bios_update_obj = { "info": bios.get("info", {}) }

        # ... do stuff ...

        if not byc["test_mode"]:
            #bios_coll.update_one( { "_id": _id }, { '$set': bios_update_obj }  )
            continue
        else:
            prjsonnice(bios_update_obj)

        ####################################################################
        ####################################################################
        ####################################################################

    bar.finish()

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
