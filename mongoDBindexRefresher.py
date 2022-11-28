#!/usr/local/bin/python3

import re, json, yaml, statistics, sys, datetime, time
from os import path, environ, pardir
from pymongo import MongoClient, GEOSPHERE
from random import sample as randomSamples
from progress.bar import Bar

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""
"""

################################################################################
################################################################################
################################################################################

def main():

	mongodb_index_refresher()

################################################################################

def mongodb_index_refresher():

	initialize_service(byc, "biosamples_refresher")
    
	dt_m = byc["datatable_mappings"]
	b_rt_s = byc["beacon_mappings"]["response_types"]

	ds_id = "progenetix"
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

	"""

#!/bin/bash

# Script to ensure MongoDB indexing
# Please modify for additional collections and changing attributes.

#
# progenetix publications
#

for db in progenetix
do
    for dbcoll in publications
    do
        echo "=> index for $db.$dbcoll.provenance.geo_location.geometry"
        mongo $db --eval "db.$dbcoll.createIndex( { \"provenance.geo_location.geometry\" : \"2dsphere\" } )"
    done
done


for field in \"id\" \"provenance.geo.city\" \"counts.ccgh\" \"counts.acgh\" \"counts.wes\" \"counts.wgs\" \"counts.ngs\" \"counts.genomes\" \"external_references.type.id\"
do
  mongo progenetix --eval "db.publications.createIndex( { $field : 1 } )"
done

#
# progenetix genes
#
for field in symbol accession_version reference_name start end gene_id swiss_prot_accessions
do
  mongo progenetix --eval "db.genes.createIndex( { $field : 1 } )"
done


for db in progenetix
do
	for field in updated \"variant_state.id\" variant_type reference_bases callset_id alternate_bases individual_id biosample_id \"type\" \"location.sequence_id\" \"location.interval.start.value\" \"location.interval.end.value\" \"info.var_length\" variant_internal_id
	do
		echo "=> index for $db.variants.$field"
		mongo $db --eval "db.variants.createIndex( { $field : 1 } )"
	done

    for dbcoll in biosamples callsets individuals
	do
		echo "=> index for $db.$dbcoll.id"
		mongo $db --eval "db.$dbcoll.createIndex( { 'id' : 1 }, { 'unique': true } )"
	
		for field in id \"external_references.id\" \"external_references.label\" description \"provenance.geo_location.properties.city\" \"provenance.geo_location.properties.country\" individual_id age_at_collection biosample_status.id
		do
			echo "=> index for $db.biosamples.$field"
			mongo $db --eval "db.biosamples.createIndex( { $field : 1 } )"
		done
	
	done

	for field in biosample_id individual_id
	do
		echo "=> index for $db.callsets.$field"
		mongo $db --eval "db.callsets.createIndex( { $field : 1 } )"
	done

	
    for dbcoll in collations
	do
		for field in child_terms id label count collation_type namespace_prefix
		do
			echo "=> index for $db.collations.$field"
			mongo $db --eval "db.collations.createIndex( { $field : 1 } )"
		done	
	done

	for dbcoll in biosamples individuals callsets
	do
		echo "=> index for $db.$dbcoll.provenance.geo_location.geometry"
		mongo $db --eval "db.$dbcoll.createIndex( { \"provenance.geo_location.geometry\" : \"2dsphere\" } )"
	done
done

	
for db in progenetix 
do
	for dbcoll in biosamples
	do
		echo "=> index for $db.$dbcoll histologies etc."
		for field in \"cohorts.id\" \"histological_diagnosis.id\" \"sampled_tissue.id\" \"icdo_morphology.id\" \"icdo_topography.id\" \"pathological_tnm_findings.id\" \"tumor_grade.id\" \"pathological_stage.id\"
		do
			echo "=> index for $db.$dbcoll.$field"
			mongo $db --eval "db.$dbcoll.createIndex( { $field : 1 } )"
		done
	done
done

for db in progenetix
do
	for dbcoll in individuals
	do
		echo "=> index for $db.$dbcoll histologies etc."
		for field in \"index_disease.disease_code.id\" \"auxiliary_disease.disease_code.id\" 
		do
			echo "=> index for $db.$dbcoll.$field"
			mongo $db --eval "db.$dbcoll.createIndex( { $field : 1 } )"
		done
	done
done


for db in progenetix
do
    for dbcoll in geolocs
    do
        echo "=> index for $db.$dbcoll.geo_location.geometry"
        mongo $db --eval "db.$dbcoll.createIndex( { \"geo_location.geometry\" : \"2dsphere\" } )"
    done
done

	"""

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
