#!/usr/local/bin/python3
from pymongo import MongoClient
from os import path, environ, pardir
import sys, datetime, time
from progress.bar import Bar
from copy import deepcopy as deepcopy

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""
## Usage history

### 2022-03-04

* 2966 callsets w/o variants |################################| 100% of 2966
* ===> Inserted 520566 variants

"""

################################################################################

def main():

	callsets_variants_importer()

################################################################################

def callsets_variants_importer():

	initialize_service(byc, "callsets_modifier")
	get_args(byc)
	set_processing_modes(byc)
	select_dataset_ids(byc)
	filters_from_args(byc)
	parse_variant_parameters(byc)
	generate_genomic_intervals(byc)

	processsed_root = byc["args"].source # "~/switchdrive/baudisgroup/2022-arrayexpress-reimport"
	var_temp = object_instance_from_schema_name(byc, "pgxVariant", "properties")
	prjsonnice(var_temp)

	if not byc["args"].inputfile:
		print("No inputfile file specified => quitting ...")
		exit()
	if not byc["args"].source:
		print("No source directory specified => quitting ...")
		exit()

	v_defs = byc["variant_definitions"]

	if byc["test_mode"]:
		max_count = 10
	else:
		max_count = int(byc["args"].randno)

	csids, id_types = read_tsv_to_dictlist(byc["args"].inputfile, max_count)
	ds_id = "progenetix"

	if len(byc["filters"]) > 0:
		f_l = [ f["id"] for f in byc["filters"] ]
		print(f_l)
		csids = list(filter(lambda cs: cs["callset_id"] in f_l, csids))

	data_client = MongoClient( )
	data_db = data_client[ ds_id ]
	bios_coll = data_db[ "biosamples" ]
	ind_coll = data_db[ "individuals" ]
	cs_coll = data_db[ "callsets" ]
	v_coll = data_db[ "variants" ]

	cs_no_v = []
	v_all_no = 0

	for cs_i in csids:
		csid = cs_i["callset_id"]
		if byc["test_mode"]:
			print("{}: {}".format(csid, cs_i["series_accession__id"]))

		v = v_coll.find_one({"callset_id": csid})
		# print(v)
		if not v:
			cs_no_v.append(cs_i)
		elif byc["update_mode"]:
			cs_no_v.append(cs_i)

	cs_no_v_no = len(cs_no_v)

	bar = Bar("{} callsets will receive variants".format(cs_no_v_no), max = cs_no_v_no, suffix='%(percent)d%%'+" of "+str(cs_no_v_no) )

	for cs_i in cs_no_v:

		csid = cs_i["callset_id"]

		bar.next()

		cs = cs_coll.find_one({"id": csid})
		if not cs:
			print("...no callset for {}".format(csid))
			continue

		ser_id = re.sub("arrayexpress:", "", cs_i["series_accession__id"])

		vs_path = path.join(processsed_root, ser_id, cs_i["_progenetix_experiment_id"], "variants.json")


		try:
			with open( vs_path, "r") as vip:
				vs_d = json.loads( vip.read() )
		except (FileNotFoundError, IOError) as e:
			print(e)
			continue

		if byc["update_mode"]:
			if byc["test_mode"]:
				v_d_no = len( v_coll.distinct("_id", {"callset_id": csid}) )
				print("\n! This would delete {} variants for {}; new {}".format(v_d_no, csid, len(vs_d)))
			else:
				v_coll.delete_many({"callset_id": csid})

		for v_in in vs_d:
			v = deepcopy(var_temp)
			v.update({
				"variant_internal_id": v_in["digest"],
				"callset_id": cs["id"],
				"start": int(v_in.get("start_min", None)),
				"end": int(v_in.get("end_max", None)),
				"variant_state": v_defs["cnv_iscn_defs"][ v_in["variant_type"] ]["variant_state"],
				"updated": datetime.datetime.now().isoformat()
				})
			for k in ["variant_type", "reference_name"]:
				v.update({k: v_in[k]})
			for k in ["biosample_id"]:
				v.update({k: cs[k]})
			v["info"].update({
				"var_length": v["end"]-v["start"],
				"cnv_value": v_in["info"].get("value", None),
				"probe_number": int(v_in["info"].get("probes", None))
				})

			if byc["test_mode"]:
				prjsonnice(v)
			else:
				v_u = v_coll.insert_one(v)
				v_id = "pgxvar-{}".format(v_u.inserted_id)
				v_coll.update_one( {"_id": v_u.inserted_id}, {"$set": { "id": v_id } })
				# prjsonnice(v_id)

				v_all_no += 1

		# updating the callset maps
		print("\n"+csid)
		cs_update_obj = { "info": cs["info"] }
		maps, cs_cnv_stats, cs_chro_stats = interval_cnv_arrays(v_coll, { "callset_id": csid }, byc)
		cs_update_obj.update({"cnv_statusmaps": maps})
		cs_update_obj.update({"cnv_stats": cs_cnv_stats})
		cs_update_obj.update({"cnv_chro_stats": cs_chro_stats})
		cs_update_obj.update({ "updated": datetime.datetime.now().isoformat() })

		if not byc["test_mode"]:
			cs_coll.update_one( { "_id": cs["_id"] }, { '$set': cs_update_obj }  )
		else:
			print(json.dumps(camelize(maps), sort_keys=True, default=str))

	bar.finish()
	print("===> Inserted {} variants".format(v_all_no))

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
	main()
