#!/usr/local/bin/python3
from pymongo import MongoClient
from os import path, environ, pardir
import sys, datetime
from progress.bar import Bar

import pandas as pd

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""
## Usage history

### 2022-02-17 Import of array metadata from 2 arrayExpress imports

* `./callsetsArrayexpressSampleMapper.py -i ~/groupbox/arraymapMetadata/2020-09-02-arrayExpress/arrayexprs_fixed_cutcol.tsv`
* `./callsetsArrayexpressSampleMapper.py -i ~/groupbox/arraymapMetadata/2020-10-12-arrayExpress/new_arrayexpress_samples_cancer_2nd_batch_uniq.tsv`

This resulted in the platform and accession labeling of 1074 callsets'

"""

################################################################################
def main():

	callsets_modifier()

################################################################################

def callsets_modifier():

	initialize_bycon_service(byc)
	get_args(byc)
	set_processing_modes(byc)

	select_dataset_ids(byc)

	if not byc["args"].inputfile:
		print("No inputfile file specified => quitting ...")
		exit()

	metapath = byc["args"].inputfile
	platform_rename = {
		'GPL2005': 'Xba',
		'GPL2004': 'Hind',
		'GPL3720': 'Sty',
		'GPL3718': 'Nsp',
		'GPL16131': 'CytoScanHD',
		'GPL18637': 'CytoScan750K',
		'GPL6801': 'SNP6',
		'GPL6894': 'SNP5',
		'GPL2641': 'Xba142'
	}

	pre = "GenomePlatforms"
	ds_id = "progenetix"
	f_d = byc["filter_definitions"][pre]

	pre_h_f = path.join( parent_path, "byconeer", "rsrc", pre, "numbered_hierarchies.tsv" )
	hier = hierarchy_from_file(ds_id, pre, pre_h_f, byc)
	c_l_s = {}
	for c_id, c_v in hier.items():
		r = "{}{}".format( f_d["reference"]["root"], c_id)
		c_l_s.update({ c_id: {
			"description": c_v["label"],
			"label": c_v["label"],
			"reference": r
			}
		})

	data_client = MongoClient( )
	data_db = data_client[ ds_id ]
	cs_coll = data_db[ "callsets" ]
	bios_coll = data_db[ "biosamples" ]

	### read in meta table
	mytable = pd.read_csv(metapath, sep = '\t', dtype = {'AGE': str, 'YEAR': str, 'PMID': str})
	mytable = mytable.where((pd.notnull(mytable)), None) ## convert pd.nan to None
	no_row = mytable.shape[0]

	bar = Bar("{} experiments".format(no_row), max = no_row, suffix='%(percent)d%%'+" of "+str(no_row) )

	for row in mytable.itertuples():

		if not byc["test_mode"]:
			bar.next()

		if hasattr(row, "status"):
			status = row.status
			if status.startswith('excluded'):
				continue
		ser = row.SERIESID
		arr = row.UID.replace("-", "_").replace(' ', '_').replace('.CEL','')
		plf = row.PLATFORMID
		apl = row.PLATFORM
		pf_apl = "arrayexpress:"+apl
		ae_s_id = "arrayexpress:{}-{}".format(ser, arr)
		callset_id_leg = "pgxcs::{}::AE-{}-{}-{}".format(ser, ser.replace('E-','').replace('-','_'), arr, platform_rename[plf])
		# print(callset_id_leg, row.PLATFORM)
		update_obj = {
			"series_accession": { "id": "arrayexpress:"+ser, "label": "" },
			"experiment_accession": { "id": ae_s_id, "label": row.UID },
			"platform_model": { "id": pf_apl, "label": "" }
		}

		if pf_apl in c_l_s.keys():
			update_obj["platform_model"].update({"label": c_l_s[pf_apl]["label"]})

		# callset update

		cs = cs_coll.find_one({"info.legacy_id": callset_id_leg})
		if not cs:
			continue

		if not byc["test_mode"]:
			cs_coll.update_one({"_id": cs["_id"]}, {"$set": update_obj })

		bios = bios_coll.find_one({"id": cs["biosample_id"]})
		if not bios:
			continue

		# biosample update

		e_r_s = [ ]
		for e_r in bios["external_references"]:
			if not "arrayexpress" in e_r["id"]:
				e_r_s.append(e_r)
		e_r_s.append( {"id": "arrayexpress:"+ser, "label": "" } )
		e_r_s.append( {"id": ae_s_id, "label": row.UID } )
		update_obj = { "external_references": e_r_s }
		if not byc["test_mode"]:
			bios_coll.update_one({"_id": bios["_id"]}, {"$set": update_obj })

	bar.finish()

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
	main()
