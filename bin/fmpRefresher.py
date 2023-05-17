#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
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

    initialize_bycon_service(byc, "biosamples_refresher")
    get_args(byc)
    set_processing_modes(byc)

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
    parse_variant_parameters(byc)
    initialize_beacon_queries(byc)
    generate_genomic_intervals(byc)

    io_params = byc["datatable_mappings"]["io_params"][ "biosample" ]
    io_prefixes = byc["datatable_mappings"]["io_prefixes"][ "biosample" ]

    with open( path.join(byconeer_path, *byc["this_config"]["biosample_id_path"]) ) as lif:
        bs_legacy_ids = yaml.load( lif, Loader=yaml.FullLoader)

    legacy_ids = {}
    for bs_id, l_id in bs_legacy_ids.items():
        legacy_ids.update({l_id: bs_id})

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    bios_coll = data_db[ "biosamples" ]
    ind_coll = data_db[ "individuals" ]
    cs_coll = data_db[ "callsets" ]
    v_coll = data_db[ "variants" ]
    no_match = 0
    counter = 0
    max_count = 1000000

    logfile = path.splitext(byc["args"].inputfile)[0]
    if byc["test_mode"]:
        logfile += "_test.tsv"
        max_count = 2000
    else:
        logfile += "_processed.tsv"
    log = open(logfile, "w")

    fmp_samples, fmp_pmids, fieldnames = _read_samplefile(counter, max_count, byc)
    fmp_no = len(fmp_samples)
    print("=> The samplefile contains {} samples from {} studies".format(fmp_no, len(fmp_pmids)))

    missing_ids = []
    incomplete_pubs = set()

    # for f_s in fmp_samples:
    #     bios_id = legacy_ids.get(f_s["pgx_fmp_sample_id"], None)
    #     if bios_id is None:
    #         missing_ids.append(f_s["pgx_fmp_sample_id"])
    #         incomplete_pubs.add(f_s["external_references__id___PMID"])
    #         print(f_s["pgx_fmp_sample_id"])

    csv_writer = csv.DictWriter(log, delimiter="\t", fieldnames=fieldnames)
    csv_writer.writeheader()

    sel_pmids = []
    sel_fmp_samples = []

    for fmp_pmid in fmp_pmids:

        this_study = False

        f_p_samples = list(filter(lambda p: p["external_references__id___PMID"] == fmp_pmid, fmp_samples))

        for f_s in f_p_samples:
            # TODO: move away from hard coded criterium ...
            if "ISCN" in f_s["data_provenance"]:
                this_study = True
 
        if this_study is True:
            sel_pmids.append(fmp_pmid)
            for f_s in f_p_samples:
                sel_fmp_samples.append(f_s)

    sel_no = len(sel_fmp_samples)

    print("=> {} of {} studies were selected for processing".format(len(sel_pmids), len(fmp_pmids)))

    pgx_samples, pgx_pmids = _get_pmids_samples(sel_pmids, bios_coll, byc)
    print("=> The database contains {} samples for {} of those studies".format(len(pgx_samples), len(pgx_pmids)))

    for fmp_pmid in sel_pmids:

        f_p_samples = list(filter(lambda p: p["external_references__id___PMID"] == fmp_pmid, sel_fmp_samples))
        f_p_no = len(f_p_samples)

        p_p_samples = list(filter(lambda p: p["external_references__id___PMID"] == fmp_pmid, pgx_samples))
        p_p_no = len(p_p_samples)

        if p_p_no >= f_p_no:
            continue

        print("{} has {} samples in import and {} in pgx".format(fmp_pmid, f_p_no, p_p_no))

        for f_s in f_p_samples:

            legacy_id_construct = "PGX_AM_BS_"+f_s["pgx_fmp_sample_id"]
            legacy_id = legacy_ids.get(f_s["pgx_fmp_sample_id"], "___nothing___")

            pgx_s_ids = []
            for p_s in p_p_samples:
                pgx_s_ids.append(p_s["id"])

            if legacy_id in pgx_s_ids:
                continue

            pgx_s_l_ids = []
            for p_s in p_p_samples:
                pgx_s_l_ids += p_s["info"]["legacy_ids"]

            if legacy_id_construct in pgx_s_l_ids:
                continue

            no_match += 1

            # now updating sample data ...

            h_query = {
                "icdo_morphology.id": f_s["icdo_morphology__id"],
                "icdo_topography.id": f_s["icdo_topography__id"]
            }

            histomatch = bios_coll.find_one(h_query)
            if histomatch:
                try:
                    f_s.update({
                        "histological_diagnosis__id": histomatch["histological_diagnosis"].get("id", ""),
                        "histological_diagnosis__label": histomatch["histological_diagnosis"].get("label", ""),
                        "sampled_tissue__id": histomatch["sampled_tissue"].get("id", ""),
                        "sampled_tissue__label": histomatch["sampled_tissue"].get("label", "")
                    })
                except:
                    pass

            csv_writer.writerow(f_s)

        for p_s in p_p_samples:
            match_status = p_s.get("_progenetix_status", None)
            if not match_status:
                p_e_s = {
                    "id": p_s["id"],
                    "external_references__id___PMID": fmp_pmid,
                    "info__legacy_ids": ""
                }

                try:
                    p_e_s.update({"info__legacy_ids": ",".join(p_s["info"].get("legacy_ids", []))})
                except:
                    pass

                # TODO: Obviously should be generalized/moved/subbed
                for p, k in io_params.items():
                    if p in fieldnames:
                        v = get_nested_value(p_s, k["db_key"])
                        if isinstance(v, list):
                            p_e_s.update({ p: "::".join(v) } )
                        else:
                            if "integer" in k["type"]:
                                 p_e_s.update({ p: int(v) })
                            p_e_s.update({ p: str(v) })

                        # if len(p_e_s[p]) > 0:
                        #     print("found some {}: {}".format(p, p_e_s[p]))

                for par, d in io_prefixes.items():

                    exts = ["id", "label"]

                    if "string" in d["type"]:
                        if par in p_s:
                            for e in exts:
                                p = par+"__"+e
                                if p in fieldnames:
                                    v = p_s[par].get(e, "")
                                    p_e_s.update({ p: v })
                                    # print("found some {}: {}".format(p, p_e_s[p]))
                    elif "array" in d["type"]:
                        for pre in d["pres"]:
                            for o in p_s[par]:
                                if pre in o["id"]:
                                    for e in exts:
                                        p = par+"__"+e+"___"+pre
                                        if p in fieldnames:
                                            v = o.get(e, "")
                                            p_e_s.update({ p: v })
                                            # print("found some {}: {}".format(p, p_e_s[p]))

                csv_writer.writerow(p_e_s)

    print("=> {} / {} samples were not found".format(no_match, len(sel_fmp_samples)))
    exit()

################################################################################

def _read_samplefile(counter, max_count, byc):

    fmp_samples = []
    fmp_pmids = set()
    fieldnames = [ "id" , "info__legacy_ids" ]

    with open(byc["args"].inputfile, newline='') as csvfile:
        
        fmp_in = csv.DictReader(csvfile, delimiter="\t", quotechar='"')

        fieldnames += fmp_in.fieldnames

        for fmp_s in fmp_in:

            fmp_s = dict(fmp_s)

            if counter > max_count:
                return fmp_samples, fmp_pmids, fieldnames
            counter += 1

            pmid = fmp_s["external_references__id___PMID"]
            if not "PMID:" in pmid:
                pmid = "PMID:"+pmid

            fmp_s.update({"external_references__id___PMID":pmid})

            fmp_pmids.add(pmid)
            fmp_samples.append(dict(fmp_s))

    fieldnames += [
        "histological_diagnosis__id",
        "histological_diagnosis__label",
        "sampled_tissue__id",
        "sampled_tissue__label",
        "_progenetix_status"
    ]

    return fmp_samples, list(fmp_pmids), fieldnames

################################################################################

def _get_pmids_samples(fmp_pmids, bios_coll, byc):

    pgx_samples = []
    pgx_pmids = set()

    for pmid in fmp_pmids:
        for bios in bios_coll.find( {"external_references.id":pmid }):
            pgx_pmids.add(pmid)
            bios.update({"external_references__id___PMID":pmid})
            pgx_samples.append(bios)

    return pgx_samples, pgx_pmids

################################################################################

def _ind_update_survival(ind, fmp_s, ind_coll, byc):

    print("\n{} => {}".format(fmp_s["death"], ind["vital_status"]["status"]))

    ind_up = {}

    # old data for modding
    ind_up.update({
        "vital_status": ind.get( "vital_status", { "status":"UNKNOWN_STATUS" } ),
        "diseases": ind.get("diseases", [{}])
        })

    if fmp_s["followup_state__id"] == "EFO:0030049":
        ind_up.update({"vital_status":{"status":"DECEASED"}})
        ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030049", "label":"dead (follow-up status)"}})

    if fmp_s["followup_state__id"] == "EFO:0030050":
        if not "DECEASED" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"DECEASED"}})
        ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030050", "label":"death from disease"}})

    if fmp_s["followup_state__id"] == "EFO:0030041":
        if not "ALIVE" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"ALIVE"}})
            ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030041", "label":"alive (follow-up status)"}})

    if fmp_s["followup_state__id"] == "EFO:0030048":
        if not "ALIVE" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"ALIVE"}})
        ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030048", "label":"alive in complete remission"}})

    if fmp_s["followup_state__id"] == "EFO:0030042":
        if not "ALIVE" in ind["vital_status"]["status"]:
            ind_up.update({"vital_status":{"status":"ALIVE"}})
        ind_up["diseases"][0].update({"followup_state":{"id":"EFO:0030042", "label":"alive with disease"}})

    if not byc["test_mode"]:
        ind_coll.update_one({"_id": ind["_id"]}, {"$set": ind_up })

################################################################################

def _ind_update_sex(ind, fmp_s, ind_coll, byc):

    ind_up = { "sex": ind.get("sex", { "id": 'PATO:0020000', "label": 'genotypic sex' }) }

    if "PATO:0020002" in fmp_s["sex_id"]:
        ind_up["sex"].update({"id": 'PATO:0020002', "label": 'female genotypic sex'})
    elif "PATO:0020001" in fmp_s["sex_id"]:
        ind_up["sex"].update({"id": 'PATO:0020001', "label": 'male genotypic sex'})

    if not byc["test_mode"]:
        ind_coll.update_one({"_id": ind["_id"]}, {"$set": ind_up })

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
