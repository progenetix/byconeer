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
## `variantsExporter`

* variantsExporter.py -f "icdot-C71.6" -d progenetix -o ~/Downloads/test.pgxseg
"""

################################################################################

def main():

    variants_checker()

################################################################################

def variants_checker():

    initialize_service(byc, "variants_exporter")
    get_args(byc)
    set_test_mode(byc)

    select_dataset_ids(byc)
    check_dataset_ids(byc)

    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d ...")
        exit()
    if len(byc["dataset_ids"]) > 1:
        print("Only one dataset can be provided with -d ...")
        exit()

    ds_id = byc[ "dataset_ids" ][ 0 ]

    mongo_client = MongoClient( )
    vs_coll = mongo_client[ ds_id ][ "variants" ]
    bs_coll = mongo_client[ ds_id ][ "biosamples" ]
    cs_coll = mongo_client[ ds_id ][ "callsets" ]

    bs_ids = vs_coll.distinct("biosample_id")
    cs_ids = vs_coll.distinct("callset_id")

    all_vs_no = vs_coll.estimated_document_count()
    all_bs_no = len(bs_ids)
    all_cs_no = len(cs_ids)

    ############################################################################

    bar = Bar("Testing {} variants for {} biosamples".format(all_vs_no, all_bs_no), max = all_bs_no, suffix='%(percent)d%%' )

    bs_ids_missing = []
    vs_missing_bs = []

    all_v_no = 0
    for bs_id in bs_ids:
        bs = bs_coll.find_one({ "id": bs_id })
        if not bs:
            bs_ids_missing.append(bs_id)
        bar.next()
    bar.finish()

    missing_bs_no = len(bs_ids_missing)

    print("=> {} biosamples were missing for variants...".format(missing_bs_no))

    print("\n"+"\n".join(bs_ids_missing)+"\n")

    for bs_id in bs_ids_missing:
        for v in vs_coll.find({"biosample_id": bs_id}):
            vs_missing_bs.append(v)

    missing_bs_v_no = len(vs_missing_bs)

    print("=> {} variants belonged to those biosamples...".format(missing_bs_v_no))

    ############################################################################

    bar = Bar("Testing {} variants for {} callsets".format(all_vs_no, all_cs_no), max = all_cs_no, suffix='%(percent)d%%' )

    cs_ids_missing = []
    vs_missing_cs = []

    all_v_no = 0
    for cs_id in cs_ids:
        cs = cs_coll.find_one({ "id": cs_id })
        if not cs:

            v = vs_coll.find_one({"callset_id":cs_id})

            v_bs_id = v.get("biosample_id", "__nix__")
            if v_bs_id not in bs_ids_missing and "__nix__" not in v_bs_id:

                print("\n"+'¯\_(ツ)_/¯ creating cs for {}'.format(v_bs_id))






            else:



                cs_ids_missing.append(cs_id)




        bar.next()
    bar.finish()

    missing_cs_no = len(cs_ids_missing)

    print("=> {} callsets were missing for variants...".format(missing_cs_no))

    # print("\n"+"\n".join(cs_ids_missing)+"\n")

    for cs_id in cs_ids_missing:
        for v in vs_coll.find({"callset_id": cs_id}):
            vs_missing_cs.append(v)

        vs_coll.delete_many({"callset_id":cs_id})

    missing_cs_v_no = len(vs_missing_cs)

    print("=> {} variants belonged to those callsets were deleted".format(missing_cs_v_no))




    exit()

    # f = open(byc["args"].outfile, "w")

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
