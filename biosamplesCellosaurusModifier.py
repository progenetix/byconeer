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

    callsets_modifier()

################################################################################

def callsets_modifier():

    initialize_service(byc)
    get_args(byc)
    set_test_mode(byc)

    select_dataset_ids(byc)
    check_dataset_ids(byc)

    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d ...")
        exit()

    pre = "cellosaurus"

    f_d = byc["filter_definitions"][pre]

    bios_q = { "external_references.id":{"$regex": pre} }
    coll_q = { "id":{"$regex": pre} }

    byc["config"]["dataset_ids"] = [ "progenetix" ]

    print(bios_q)

	############################################################################

    for ds_id in byc["config"]["dataset_ids"]:

	    data_client = MongoClient( )
	    data_db = data_client[ ds_id ]
	    cs_coll = data_db[ "callsets" ]
	    bios_coll = data_db[ "biosamples" ]
	    coll_coll = data_db[ "collations" ]

	    pre_h_f = path.join( parent_path, "byconeer", "rsrc", pre, "numbered-hierarchies.tsv" )
	    hier = hierarchy_from_file(ds_id, pre, pre_h_f, byc)
	    c_l_s = {}
	    for c_id, c_v in hier.items():
	    	r = "{}{}".format( f_d["url"], re.sub(pre+':', "", c_id) )
	    	c_l_s.update({ c_id: {
	    		"description": c_v["label"],
	    		"reference": r
	    		}
	    	})

	    # coll_labels = coll_coll.find(coll_q, { "_id": 0, "id": 1, "label":1, "reference": 1})
	    # cs_no = cs_coll.estimated_document_count()

	    no = len(list(bios_coll.find(bios_q)))
	    bios_cursor = bios_coll.find(bios_q)

	    bar = Bar("{} matched biosamples from {}".format(no, ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )

	    for bios in bios_cursor:

	    	bios_id = bios["id"]

	    	if not byc["test_mode"]:
	    		bar.next()

	    	# cs = cs_coll.find_one({"biosample_id": bios_id})
	    	# if not cs:
	    	# 	continue

	    	

	    	e_r_s = []

	    	for e_r in bios["external_references"]:
	    		c_id = e_r["id"]
	    		if "cellosaurus" in c_id:
	    			if "unknown" in c_id:
	    				continue
	    			c_id = re.sub("cvcl", "CVCL", c_id)
	    			if c_id in c_l_s.keys():
	    				c_d = c_l_s[c_id]
	    				for k, v in c_d.items():
	    					e_r.update({k: v})
	    				e_r.pop('label', None)

	    			e_r_s.append(e_r)
	    		else:
	    			e_r_s.append(e_r)

	    	update_obj = { "external_references": e_r_s }

	    		# if "geo:GPL" in e_r["id"]:
	    		# 	p_l = ""
	    		# 	if e_r["id"] in c_l_s.keys():
	    		# 		p_l = c_l_s[ e_r["id"] ]
	    		# 	update_obj.update({"platform_model": {"id": e_r["id"], "label":p_l}})
	    		# if "geo:GSM" in e_r["id"]:
	    		# 	update_obj.update({"experiment_accession": {"id": e_r["id"], "label":e_r.get("label", "")}})
	    		# if "geo:GSE" in e_r["id"]:
	    		# 	update_obj.update({"series_accession": {"id": e_r["id"], "label":e_r.get("label", "")}})
	    		# if "label" in e_r:
		    	# 	if "TCGA sample_id" in e_r["label"]:
		    	# 		update_obj.update({"experiment_accession": {"id": e_r["id"], "label":e_r.get("label", "")}})
		    	# 		update_obj.update({"platform_model": {"id": "geo:GPL6801", "label":"[GenomeWideSNP_6] Affymetrix Genome-Wide Human SNP 6.0 Array"}})
		    	# 	if "TCGA" in e_r["label"] and "project" in e_r["label"]:
		    	# 		update_obj.update({"series_accession": {"id": e_r["id"], "label":e_r.get("label", "")}})

	    	if not byc["test_mode"]:
    			bios_coll.update_one({"_id": bios["_id"]}, {"$set": update_obj })


	    bar.finish()

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
