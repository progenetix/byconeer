#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from pymongo import MongoClient
import statistics
from progress.bar import Bar
import csv

dir_path = path.dirname( path.abspath(__file__) )
byconeer_path = path.join( dir_path, pardir )
parent_path = path.join( byconeer_path, pardir )
sys.path.append( parent_path )

from byconeer import *
from bycon import *

################################################################################
################################################################################
################################################################################

def main():

    publication_tagger()

################################################################################

def publication_tagger():

    initialize_bycon_service(byc)
    run_beacon_init_stack(byc)

    if len(byc["dataset_ids"]) > 1:
        print("Please give only one dataset using -d")
        exit()

    byc.update({
        "queries": {
            "biosamples": {
                "external_references.id": re.compile( r'cellosaurus:', re.IGNORECASE )
            }
        }
    })

    ds_id = byc["dataset_ids"][0]

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    bios_coll = data_db[ "biosamples" ]
    pub_coll = data_client["progenetix"][ "publications" ]

    set_processing_modes(byc)
    # execute_bycon_queries( ds_id, byc )

    bs_no = bios_coll.count_documents(byc["queries"]["biosamples"])

    bar = Bar("Found {} cell line samples".format(bs_no), max = bs_no, suffix="%(percent)d%%"+" of "+str(bs_no) ) 

    pubs = {}

    for bios in bios_coll.find(byc["queries"]["biosamples"]):

        _id = bios["_id"]

        bar.next()
        e_r_s = bios.get("external_references", [])
        pmid = False
        for e_r in e_r_s:
            if "PMID" in e_r["id"]:
                pmid = e_r["id"]
                if pmid in pubs.keys():
                    continue
                else:
                    pubs.update({pmid:{"cell_lines":{}}})
                    continue
        if not pmid:
            continue
        for e_r in e_r_s:
            if "cellosaurus" in e_r["id"]:
                if e_r["id"] in pubs[pmid]["cell_lines"].keys():
                    continue
                pubs[pmid]["cell_lines"].update({
                    e_r["id"]: e_r
                })

    for p_k, p_v in pubs.items():
        if byc["update_mode"] is True:
            p_e = pub_coll.update_one(
                {"id": p_k},
                {"$set": {
                        "cell_lines": list(p_v["cell_lines"].values()),
                        "has_cell_lines": True
                    }
                }
            )
            # print(p_e)

        else:
            p_e = pub_coll.find_one({"id": p_k})
            if p_e is False:
                print(p_k)
            prjsonnice(list(p_v["cell_lines"].values()))

    bar.finish()

    print("Found {} publications for {} cell lines...".format(len(pubs), bs_no))

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
