#!/usr/local/bin/python3
from pymongo import MongoClient
from os import path, environ, pardir
import re, sys, datetime
from progress.bar import Bar

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

################################################################################
def main():

    biosamples_refresher()

################################################################################

def biosamples_refresher():

    initialize_bycon_service(byc)
    
    ds_id = "cellz"

	############################################################################

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    cs_coll = data_db[ "callsets" ]
    bios_coll = data_db[ "biosamples" ]
    coll_coll = data_db[ "collations" ]

    bios_q = {}

    no = len(list(bios_coll.find(bios_q)))
    bios_cursor = bios_coll.find(bios_q)

    bar = Bar("{} matched biosamples from {}".format(no, ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )

    for bios in bios_cursor:

    	bios_id = bios["id"]

    	if not byc["test_mode"]:
    		bar.next()

    	cl_obj = {
    		"id": "NA",
    		"label": "NA",
    		"reference": "",
    		"synonyms": [],
    		"notes": ""
    	} 

    	for e_r in bios["external_references"]:
    		c_id = e_r["id"]
    		if "cellosaurus" in c_id:
    			if "unknown" in c_id:
    				continue
    			cl_obj.update({
    				"id": c_id,
    				"label": e_r.get("description", ""),
    				"reference": e_r.get("reference", ""),
    				"synonyms": [c_id, re.sub("cellosaurus:", "", c_id )]
    				})
    			if len(cl_obj["label"]) > 0:
    				cl_obj["synonyms"].append(cl_obj["label"])
    			continue

    	if not byc["test_mode"]:
    		bios_coll.update_one({"_id": bios["_id"]}, {"$set": {"cellline_info": cl_obj} })

    bar.finish()

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
