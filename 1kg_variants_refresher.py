#!/usr/local/bin/python3

import re, json, yaml, sys, datetime, statistics
from os import path, environ, pardir
from pymongo import MongoClient
from isodate import date_isoformat
from progress.bar import Bar

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""
## `biosamplesRefresher`

"""

################################################################################
################################################################################
################################################################################

def main():

    variants_refresher()

################################################################################

def variants_refresher():

    initialize_service(byc)
    select_dataset_ids(byc)
    print(byc)
    
    if len(byc["dataset_ids"]) < 1:
        print("No existing dataset was provided with -d ...")
        exit()

    if byc["test_mode"]:
        max_count = 10
    else:
        max_count = int(byc["args"].randno)

    mongo_client = MongoClient( )

    # for ds_id in byc["dataset_ids"]:

    #     data_db = mongo_client[ ds_id ]
    #     var_coll = data_db[ "variants" ]
    #     bios_coll = data_db[ "biosamples" ]
    #     no =  bios_coll.estimated_document_count()
    #     if max_count < 1:
    #         max_count = no
    #     count = 0
    #     bar = Bar("{} vars".format(ds_id), max = no, suffix='%(percent)d%%'+" of {} biosamples".format(no) )
    #     for bios in bios_coll.find({}):
    #         bios_id = bios["id"]
    #         ind_id = bios.get("individual_id", "")

    #         for v_p in var_coll.find({"biosample_id": bios_id}):
    #             var_coll.update_one( { "_id": v_p["_id"] }, { '$set': {"individual_id": ind_id} }  )            
    #         bar.next()
    #     bar.finish()

    # exit()

    v_d = byc["variant_definitions"]

    for ds_id in byc["dataset_ids"]:

        v_short = 0
        v_no_type = 0

        data_db = mongo_client[ ds_id ]
        var_coll = data_db[ "variants" ]
        bios_coll = data_db[ "biosamples" ]
        no =  var_coll.estimated_document_count()
        if max_count < 1:
            max_count = no
        count = 0
        bar = Bar("{} vars".format(ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )
        for v_p in var_coll.find({"biosample_id" : "onekgbs-HG00100"}):

            if count > max_count:
                break
            count += 1

            update_obj = vrsify_variant(v_p, v_d)
            update_obj.update({
                "variant_internal_id": variant_create_digest(v_p, byc),
                "updated": datetime.datetime.now().isoformat()
            })

            update_obj.pop("digset", None) # this was a buggy remnant

            # for p in ["callset_id", "biosample_id", "variant_internal_id", "variant_type", "variant_state", "reference_bases", "alternate_bases", "info"]:
            #     p_v = v_p.get(p, None)
            #     if p_v is not None:
            #         update_obj.update({p:p_v})

            if not "start" in v_p:
                break

            if not byc["test_mode"]:
                var_coll.update_one( { "_id": v_p["_id"] }, { '$set': update_obj }  )
            else:
                prjsonnice(update_obj)
            
            bar.next()
        bar.finish()

################################################################################

################################################################################

def refseq_from_chro(chro, v_d):

    chro = re.sub("chr", "", chro) # just making sure ...

    if not chro in v_d["chro_refseq_ids"]:
        return chro

    return v_d["chro_refseq_ids"][ chro ]

################################################################################

def efo_to_vrs_class(efo_id, v_d):
    
    efo_vrs = v_d["efo_vrs_map"]
    if efo_id in efo_vrs:
        return efo_vrs[ efo_id ]["relative_copy_class"]

    return None

################################################################################

def vrsify_variant(v, v_d):

    t = v.get("type", None)
    v_t = v.get("variant_type", None)
    a_b = v.get("alternate_bases", None)

    if t is None:
        if a_b is not None:
            t = "Allele"
        elif v_t in ["DUP", "AMP", "DEL", "HOMODEL"]:
            t = "RelativeCopyNumber"
    print(t)

    if "Allele" in t:
        return vrsify_snv(v, v_d)
    else:
        return vrsify_cnv(v, v_d)

################################################################################

def vrsify_cnv(v, v_d):

    v_t = v.get("type", "RelativeCopyNumber")
    ref_id = refseq_from_chro(v["reference_name"], v_d)
    start_v = int(v["start"])
    end_v = int( v.get("end",start_v + 1 ) )

    vrs_class = None
    if "variant_state" in v:
        vrs_class = efo_to_vrs_class(v["variant_state"].get("id", None), v_d)

    vrs_v = {
                "type": v_t,
                "relative_copy_class": vrs_class,
                "location":{
                    "sequence_id": ref_id,
                    "type": "SequenceLocation",
                    "interval": {
                        "start": {
                            "type": "Number",
                            "value": start_v
                        },
                        "end": {
                            "type": "Number",
                            "value": end_v
                        }
                    }
                }
            }

    return vrs_v

################################################################################

def vrsify_snv(v, v_d):

    v_t = v.get("type", "Allele")
    ref_id = refseq_from_chro(v["reference_name"], v_d)
    alt = v.get("alternate_bases", None)

    if alt is None:
        return {}

    ref_l = int(len(v.get("reference_bases", "")))
    alt_l = int(len(v.get("alternate_bases", "")))
    l = alt_l - ref_l + 1

    ref_id = refseq_from_chro(v["reference_name"], v_d)
    start_v = int(v["start"])
    end_v = int( start_v + l )

    vrs_v = {
                "type": v_t,
                "state": {
                    "type": "LiteralSequenceExpression",
                    "sequence": alt
                },
                "location":{
                    "sequence_id": ref_id,
                    "type": "SequenceLocation",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {
                            "type": "Number",
                            "value": start_v
                        },
                        "end": {
                            "type": "Number",
                            "value": end_v
                        }
                    }
                }
            }

    return vrs_v

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
