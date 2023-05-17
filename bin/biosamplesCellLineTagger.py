#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from pymongo import MongoClient
import statistics
from progress.bar import Bar
import csv


# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

################################################################################
################################################################################
################################################################################

def main():

    biosamples_tagger()

################################################################################

def biosamples_tagger():

    initialize_bycon_service(byc)
    select_dataset_ids(byc)
    parse_filters(byc)
    parse_variant_parameters(byc)
    initialize_beacon_queries(byc)
    generate_genomic_intervals(byc)
    set_processing_modes(byc)

    byc["dataset_ids"] = ['progenetix']

    if len(byc["dataset_ids"]) != 1:
        print("No single, existing dataset was provided with -d ...")
        exit()

    ds_id = byc["dataset_ids"][0]

    if not byc["args"].inputfile:
        print("No inputfile file specified => quitting ...")
        exit()

    # TODO: Move  to sub ...
    i_t = []
    f_p = byc["args"].inputfile
    with open(f_p) as f:
       rd = csv.reader(f, delimiter="\t", quotechar='"')
       for i, row in enumerate(rd):
           if i > 0:
               i_t.append(row)

    cvcl_ls = {}

    for row in i_t:
        if row[0].startswith('#'):
            continue
        cvcl_ls.update({row[0]: row[1]})

    row_no = len(cvcl_ls)
    not_found = []

    data_client = MongoClient( )
    bios_coll = data_client[ ds_id ][ "biosamples" ]
 
    bar = Bar("Reading in metadata table", max = row_no, suffix="%(percent)d%%"+" of "+str(row_no) )
    for cl_id, cl_l in cvcl_ls.items():

        bios = bios_coll.find({"external_references.id":cl_id})

        for bs in bios:

            e_r_n = []
            for e_r in bs["external_references"]:
                if cl_id in e_r["id"]:
                    e_r.update({"label":cl_l})
                e_r_n.append(e_r)
 
            if byc["update_mode"] is True:
                bios_coll.update_one(
                    {"_id":bs["_id"]},
                    {"$set": { "external_references": e_r_n} }
                )
                bar.next()
            else:
                print(e_r_n)

    bar.finish()


################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
