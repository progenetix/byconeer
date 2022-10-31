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

"""
## `variantsExporter`

* variantsExporter.py -f "icdot-C71.6" -d progenetix -o ~/Downloads/test.pgxseg
"""

################################################################################
################################################################################
################################################################################

def _get_args(byc):

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datasetIds", help="dataset")
    parser.add_argument("-f", "--filters", help="prefixed filter values, comma concatenated")
    parser.add_argument("--value-only", dest='value_only', action='store_true', help="only output variants with values")
    parser.set_defaults(value_only=False)
    byc.update({ "args": parser.parse_args() })

    return byc

################################################################################

def main():

    variants_exporter()

################################################################################

def variants_exporter():

    initialize_service(byc)
    _get_args(byc)

    select_dataset_ids(byc)
    
    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d ...")
        exit()
    if len(byc["dataset_ids"]) > 1:
        print("Only one dataset can be provided with -d ...")
        exit()

    ds_id = byc[ "dataset_ids" ][ 0 ]

    get_filter_flags(byc)  
    parse_filters(byc)
    generate_queries(byc)
    execute_bycon_queries( ds_id, byc )

    ds_results = byc["dataset_results"][ds_id]

    all_bs_no = ds_results["biosamples.id"]["target_count"]
    mongo_client = MongoClient( )
    vs_coll = mongo_client[ ds_id ][ "variants" ]
    bs_coll = mongo_client[ ds_id ][ "biosamples" ]

    bar = Bar("Testing {} {} samples".format(all_bs_no, ds_id), max = all_bs_no, suffix='%(percent)d%%' )

    used_bs_ids = set()
    all_v_no = 0
    for bs_id in ds_results["biosamples.id"]["target_values"]:
        for v in vs_coll.find({ "biosample_id": bs_id }):
            all_v_no += 1
            val = ""
            if "info" in v:
                if "cnv_value" in v["info"]:
                    if isinstance(v["info"]["cnv_value"],float):
                        val = v["info"]["cnv_value"]
            break
        if (not byc["args"].value_only) or (val != ""):
            used_bs_ids.add(bs_id)
        bar.next()
    bar.finish()

    used_bs_no = len(used_bs_ids)

    print("=> Writing header data for {} samples with variants...".format(used_bs_no))

    f = open(byc["args"].outputfile, "w")

    for bs_id in used_bs_ids:
        bs = mongo_client[ ds_id ][ "biosamples" ].find_one( { "id": bs_id } )
        h_line = "#biosample_id={}".format(bs_id)
        h_d = bs[ "histological_diagnosis" ]
        h_line = '{};group_id={};group_label={};NCIT::id={};NCIT::label={}'.format(h_line, h_d.get("id", "NA"), h_d.get("label", "NA"), h_d.get("id", "NA"), h_d.get("label", "NA"))
        f.write(h_line+"\n")

    bar = Bar("Variants from {} {} samples".format(used_bs_no, ds_id), max = used_bs_no, suffix='%(percent)d%%' )

    used_v_no = 0
    for bs_id in used_bs_ids:
        for v in vs_coll.find({ "biosample_id": bs_id }):
            val = ""
            if "info" in v:
                if "cnv_value" in v["info"]:
                    if isinstance(v["info"]["cnv_value"],float):
                        val = v["info"]["cnv_value"]
            if (not byc["args"].value_only) or (val != ""):
                used_v_no +=1
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bs_id, v["reference_name"], int(v["start"]), int(v["end"]), v["variant_type"], '{:.4f}'.format(val)) )

        bar.next()
    bar.finish()

    f.close()

    print("{} variants from {} biosamples were written to {}".format(used_v_no, len(used_bs_ids) ,byc["args"].outfile) )

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
