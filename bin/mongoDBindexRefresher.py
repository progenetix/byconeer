#!/usr/local/bin/python3

import re, json, yaml, statistics, sys, datetime, time
from os import path, environ, pardir
from pymongo import MongoClient, GEOSPHERE
from random import sample as randomSamples
from progress.bar import Bar

from bycon import *

"""
"""

################################################################################
################################################################################
################################################################################

def main():

	mongodb_index_refresher()

################################################################################

def mongodb_index_refresher():

	initialize_bycon_service(byc, "biosamples_refresher")
	select_dataset_ids(byc)
	
	dt_m = byc["datatable_mappings"]
	b_rt_s = byc["beacon_mappings"]["response_types"]

	ds_ids = byc.get("dataset_ids", [])

	if len(ds_ids) != 1:
		print("No single existing dataset was provided with -d ...")
		exit()

	ds_id = ds_ids[0]

	data_client = MongoClient( )
	data_db = data_client[ ds_id ]

	for r_t in ["biosample", "individual", "variant", "analysis"]:

		i_coll = data_db[ r_t ]
		io_params = dt_m["io_params"][ r_t ]

		for p_k, p_v in io_params.items():
			i_t = p_v.get("indexed", False)
			if i_t is False:
				continue
			k = p_v["db_key"]
			print('Creating index "{}" in {} from {}'.format(k, r_t, ds_id))
			m = i_coll.create_index(k)
			print(m)

		if i_coll.find({"provenance.geo_location.geometry": {"$exists": True}}):
			# print('geo exists in {} from {}'.format(r_t, ds_id))
			m = i_coll.create_index([("provenance.geo_location.geometry", GEOSPHERE)])
			print(m)

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
	main()
