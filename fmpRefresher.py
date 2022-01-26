#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from isodate import date_isoformat
from pymongo import MongoClient
import argparse
import statistics
from progress.bar import Bar
import base36, time

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

    fmp_refresher()

################################################################################

def fmp_refresher():

    initialize_service(byc, "biosamples_refresher")
    get_args(byc)
    set_test_mode(byc)

    if not byc["args"].inputfile:
        print("No inputfile file specified => quitting ...")
        exit()

    what = {
        "sex": "n",
        "survival": "n",
        "technique": "n"
    }

    ds_id = "progenetix"
    parse_filters(byc)
    parse_variants(byc)
    initialize_beacon_queries(byc)
    generate_genomic_intervals(byc)

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    bios_coll = data_db[ "biosamples" ]
    ind_coll = data_db[ "individuals" ]
    cs_coll = data_db[ "callsets" ]
    v_coll = data_db[ "variants" ]
    no_match = 0
    counter = 0
    max_count = 100000

    logfile = byc["args"].inputfile + ".log.tab"
    log = open(logfile, "w")
    

    with open(byc["args"].inputfile, newline='') as csvfile:
        
        fmp_in = csv.DictReader(csvfile, delimiter="\t", quotechar='"')

        fieldnames = fmp_in.fieldnames
        fmp_in = list(fmp_in)
        fmp_no = len(fmp_in)

        csv_writer = csv.DictWriter(log, delimiter="\t", fieldnames=fieldnames)
        csv_writer.writeheader()

        if not byc["test_mode"]:
            bar = Bar("{} samples will be looked up".format(fmp_no), max = fmp_no, suffix='%(percent)d%%'+" of "+str(fmp_no) )

        for fmp_s in fmp_in:

            if counter > max_count:
                return

            counter += 1

            if not byc["test_mode"]:
                max_count = fmp_no
                bar.next()

            legacy_id = "PGX_AM_BS_"+fmp_s["sample_id"]

            bios = bios_coll.find_one({"info.legacy_ids": legacy_id})
            if not bios:
                print("??? Not found: {}".format(legacy_id))
                no_match += 1
                csv_writer.writerow(fmp_s)
                continue

            bs_id = bios["id"]
            ind_id = bios["individual_id"]

            if "cCGH" in fmp_s["technique"]:

                technique = "cCGH"
           
                v_s = list(v_coll.find({"biosample_id": bs_id }))
                cs_coll.delete_one({"biosample_id": bs_id })
                v_coll.delete_many({"biosample_id": bs_id }) 
                cs_id = re.sub("bs", "cs", bs_id)

                # print("\n{} ({}): {}".format(fmp_s["sample_id"], bios["id"], cs_id, fmp_s["iscn_ccgh"]))
                variants = _deparse_rev_ish_CGH(bs_id, cs_id, technique, fmp_s["iscn_ccgh"])

                if not byc["test_mode"]:

                    for v in variants:
                        v_id = v_coll.insert_one(v).inserted_id
                        v_coll.update_one({"_id": v_id}, {"$set": {"id": str(v_id)}})
                        # print(v["biosample_id"])

                    cs_update_obj = {
                        "id": cs_id,
                        "biosample_id": bs_id,
                        "individual_id": ind_id,
                        "info": { "provenance": technique+" ISCN conversion" },
                        "platform_model": {"id": "EFO:0010937", "label":"comparative genomic hybridization (CGH)"},
                        "provenance": bios["provenance"]
                    }

                    maps, cs_cnv_stats = interval_cnv_arrays(v_coll, { "callset_id": cs_id }, byc)
                    cs_update_obj["info"].update({"statusmaps": maps})
                    cs_update_obj["info"].update({"cnvstatistics": cs_cnv_stats})

                    cs_coll.insert_one(cs_update_obj)

            ####################################################################

            ind = ind_coll.find_one({"id": ind_id})
            
            if "y" in what["survival"] and len(fmp_s["death"]) > 0:
                _ind_update_survival(ind, fmp_s, ind_coll, byc)

            if "y" in what["sex"] and len(fmp_s["sex"]) > 0:
                _ind_update_sex(ind, fmp_s, ind_coll, byc)
 
            cs_up = {}
            cs = cs_coll.find_one({"biosample_id": bs_id})

            if "y" in what["sex"] and "cCGH" in fmp_s["technique"]:
                if not "aCGH" in fmp_s["technique"]:
                    cs_up.update({"platform_model": {"id": "EFO:0010937", "label":"comparative genomic hybridization (CGH)"}})
                    if not byc["test_mode"]:
                        cs_coll.update_one({"_id": cs["_id"]}, {"$set": cs_up })
    if not byc["test_mode"]:
        bar.finish()

    print("=> {} / {} samples were not found".format(no_match, fmp_no))

################################################################################

def _ind_update_survival(ind, fmp_s, ind_coll, byc):

    print("\n{} => {}".format(fmp_s["death"], ind["vital_status"]["status"]))

    ind_up = {}

    # old data for modding
    ind_up.update({
        "vital_status": ind.get( "vital_status", { "status":"UNKNOWN_STATUS" } ),
        "diseases": ind.get("diseases", [{}])
        })

    if fmp_s["death"] == "1":
        if not "DECEASED" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"DECEASED"}})
            ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030049", "label":"dead (follow-up status)"}})

    if fmp_s["death"] == "DOD":
        if not "DECEASED" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"DECEASED"}})
        ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030050", "label":"death from disease"}})

    if fmp_s["death"] == "0":
        if not "ALIVE" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"ALIVE"}})
            ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030041", "label":"alive (follow-up status)"}})

    if fmp_s["death"] == "ACR":
        if not "ALIVE" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"ALIVE"}})
        ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030048", "label":"alive in complete remission"}})

    if fmp_s["death"] == "AWD":
        if not "ALIVE" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"ALIVE"}})
        ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030042", "label":"alive with disease"}})

    if not byc["test_mode"]:
        ind_coll.update_one({"_id": ind["_id"]}, {"$set": ind_up })

################################################################################

def _ind_update_sex(ind, fmp_s, ind_coll, byc):

    ind_up.update({
        "sex": ind.get("sex", { "id": 'PATO:0020000', "label": 'genotypic sex' })
        })

    if "female" in fmp_s["sex"]:
        if not "PATO:0020002" in ind["sex"]["id"]:
            print("\nnew female")
            ind_up["sex"].update({"id": 'PATO:0020002', "label": 'female genotypic sex'})
    elif "male" in fmp_s["sex"]:
        if not "PATO:0020001" in ind["sex"]["id"]:
            print("\nnew male")
            ind_up["sex"].update({"id": 'PATO:0020001', "label": 'male genotypic sex'})

################################################################################

def _deparse_rev_ish_CGH(bs_id, cs_id, technique, iscn):

    iscn = "".join(iscn.split())
    variants = []

    cb_pat = re.compile( byc["variant_definitions"]["parameters"]["cytoBands"]["pattern"] )

    for cnv_t, cnv_defs in byc["variant_definitions"]["cnv_iscn_defs"].items():

        revish = cnv_defs["info"]["revish_label"]

        cnv_dummy_value = byc["variant_definitions"]["cnv_dummy_values"][cnv_t]

        # print("parsing {}".format(cnv_t))

        iscn_re = re.compile(rf"^.*?{revish}\(([\w.,]+)\).*?$", re.IGNORECASE)
        if iscn_re.match(iscn):
            m = iscn_re.match(iscn).group(1)
            for i_v in re.split(",", m):

                if not cb_pat.match(i_v):
                    continue

                cytoBands, chro, start, end, error = bands_from_cytobands(i_v, byc)
                if len(error) > 0:
                    continue

                v = cnv_defs.copy()
                v.update({
                    "biosample_id": bs_id,
                    "callset_id": cs_id, 
                    "reference_name": chro,
                    "start": start,
                    "end": end,
                    "type": "CopyNumber",
                    "info": {
                        "var_length": (end - start),
                        "cnv_value": cnv_dummy_value,
                        "provenance": technique+" ISCN conversion"
                    },
                    "updated": date_isoformat(datetime.datetime.now())
                })

                v.update({"variant_internal_id": variant_create_digest(v)})

                variants.append(v)
                # print(v["variant_internal_id"])

                # if byc["test_mode"]:
                #     prjsonnice(v)

    return variants

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
