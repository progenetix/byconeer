#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir, mkdir
import sys, datetime
from pymongo import MongoClient
import argparse
import statistics
from progress.bar import Bar
import requests

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""

## `biosamplesRefresher`

* `biosamplesRefresher.py -d progenetix -m stats`

"""

################################################################################
################################################################################
################################################################################

def _get_args(byc):

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filters", help="prefixed filter value for ext. identifier")
    parser.add_argument('-s', '--scopes', help='scopes, e.g. "stage", comma-separated')
    byc.update({ "args": parser.parse_args() })

    return byc

################################################################################

def main():

    biosamples_refresher()

################################################################################

def biosamples_refresher():

    # TODO: Clean solution?

    initialize_bycon_service(byc)

    byc["dataset_ids"] = ["progenetix"]

    parse_filters(byc)
    parse_variant_parameters(byc)
    initialize_beacon_queries(byc)
    generate_genomic_intervals(byc)

    no_cs_no = 0
    no_stats_no = 0
    ds_id = byc["dataset_ids"][0]
    id_filter = "GSM"

    if byc["args"].filters:
        id_filter = byc["args"].filters

    bios_query = { "$or": [
        { "external_references.id": { "$regex": id_filter } },
        # { "info.legacy_ids": { "$regex": id_filter } },
        { "analysis_info.experiment_id": { "$regex": id_filter } }
    ] }

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    bios_coll = data_db[ "biosamples" ]
    cs_coll = data_db[ "callsets" ]

    bs_ids = []

    for bs in bios_coll.find (bios_query, {"id":1 } ):
        bs_ids.append(bs["id"])

    no =  len(bs_ids)

    bar = Bar("{} {} samples".format(no, ds_id), max = no, suffix='%(percent)d%%'+" of "+str(no) )

    series_re = re.compile(r'^.*?(GSE\d+?)(:?[^\d]|$)', re.IGNORECASE)
    experiment_re = re.compile(r'^.*?(GSM\d+?)(:?[^\d]|$)', re.IGNORECASE)
    platform_re = re.compile(r'^.*?(GPL\d+?)(:?[^\d]|$)', re.IGNORECASE)

    gsm_d_n = 0
    gse_d_n = 0
    gpl_d_n = 0

    for bsid in bs_ids:

        get_geosoft = True
        gsm = None
        gse = None
        gpl = None
        gsm_soft = False

        s = bios_coll.find_one({ "id": bsid })
        a_i = s.get("analysis_info", {})

        e_r_s = s.get("external_references", [])
        for e_r in e_r_s:
            if "GSM" in e_r["id"]:
                gsm = experiment_re.match(e_r["id"]).group(1)
            if "GSE" in e_r["id"]:
                gse = series_re.match(e_r["id"]).group(1)
            if "GPL" in e_r["id"]:
                gpl = platform_re.match(e_r["id"]).group(1)

        n_e_r_s = []
        for e_r in e_r_s:
            e_r_id = e_r.get("id", "")
            if not "geo:" in e_r_id:
                n_e_r_s.append(e_r)

        a_i = s.get("analysis_info", {})
        if gsm is None:
            gsm = a_i.get("experiment_id", None)
        if gse is None:
            gse = a_i.get("series_id", None)
        if gpl is None:
            gpl = a_i.get("platform_id", None)

        if gsm is None:
            print("\n!!!! No gsm file for {}".format(bsid))
            # missing_ids.append("\t".join([bsid, gse, gsm, gpl]))
            bar.next()
            continue

        gsm_soft, gsm_dl = _read_retrieve_save_gsm_soft(gsm, gse, byc)
        gsm_d_n += gsm_dl

        if not gsm_soft:
            print("\n!!!! No soft file for {}".format(gsm))
            bar.next()
            continue

        if gse is None:
            gse = _return_gse_from_gsm_soft_list(gsm_soft)
        if gpl is None:
            gpl = _return_gpl_from_gsm_soft_list(gsm_soft)

        gse_soft, gse_dl = _read_retrieve_save_gse_soft(gse, byc)
        gse_d_n += gse_dl
        gpl_soft, gpl_dl = _read_retrieve_save_gpl_soft(gpl, byc)
        gpl_d_n += gpl_dl

        gsm_e_r = _gsm_external_ref_from_gsm_soft(gsm_soft, gsm)
        gse_e_r = _gse_external_ref_from_gse_soft(gse_soft, gse)
        gpl_e_r = _gpl_external_ref_from_gse_soft(gpl_soft, gpl)

        n_e_r_s.append(gsm_e_r)
        n_e_r_s.append(gse_e_r)
        n_e_r_s.append(gpl_e_r)

        update_obj = {
            "external_references": n_e_r_s,
            "analysis_info": {
                "experiment_id": gsm_e_r.get("id", None),
                "series_id": gse_e_r.get("id", None),
                "platform_id": gpl_e_r.get("id", None),
            }
        }


        if not byc["test_mode"]:
            bios_coll.update_one({"_id": s["_id"]}, {"$set": update_obj })
        else:
            prjsoncam(update_obj)


        bar.next()

    bar.finish()

    print("Downloaded files:\n* {} GSM\n* {} GSE\n* {} GPL".format(gsm_d_n, gse_d_n, gpl_d_n))

################################################################################

def _read_retrieve_save_gsm_soft(gsm, gse, byc):

    i_d = re.sub("geo:", "", gsm)
    i_d_s = re.sub("geo:", "", gse)

    gsm_soft = None

    if gse is None:
        gsm_soft = retrieve_geosoft_file(i_d, byc)
        gse = _return_gse_from_gsm_soft_list(gsm_soft)

    gse_path = path.join( dir_path, *byc["geosoft_file_root"], "GSE", i_d_s )
    gse_file_path = path.join( gse_path, i_d_s+".txt" )
    gsm_file_path = path.join( gse_path, i_d+".txt" )

    if path.isfile(gsm_file_path):
        if gsm_soft is None:            
            gsm_soft = open(gsm_file_path).read().splitlines()
        return gsm_soft, 0
    else:
        if gsm_soft is None:
            gsm_soft = retrieve_geosoft_file(i_d, byc)
        if not path.isdir(gse_path):
            mkdir(gse_path)
        s_f = open(gsm_file_path, 'w')
        for l in gsm_soft:
            s_f.write(l)
            n = 1
        s_f.close()
        return gsm_soft, 1

    return False, 0

################################################################################

def _return_gse_from_gsm_soft_list(gsm_soft):

    gse_l = list(filter(lambda x:'Sample_series_id' in x, gsm_soft))
    gse = re.sub("!Sample_series_id = ", "", ges_l[0])
    gse = gse.strip()

    return gse

################################################################################

def _return_gpl_from_gsm_soft_list(gsm_soft):

    gpl_l = list(filter(lambda x:'Sample_platform_id' in x, gsm_soft))
    gpl = re.sub("!Sample_platform_id = ", "", gpl_l[0])
    gpl = gpl.strip()

    return gpl

################################################################################

def _gse_external_ref_from_gse_soft(gse_soft, gse):

    gse_e_r = {}

    gse_m = {
        "id": "!Series_geo_accession = ",
        "label": "!Series_title = "
    }

    for k, m in gse_m.items():
        for gse_l in gse_soft:
            if m in gse_l:
                gse_l = re.sub(m, "", gse_l)
                v = gse_l.strip()
                gse_e_r.update({k: v})
                continue

    gse_e_r.update({"reference":byc["filter_definitions"]["GEOseries"]["reference"]["root"]+gse_e_r["id"]})
    gse_e_r.update({"id":"geo:"+gse_e_r["id"]})

    return gse_e_r

################################################################################

def _gsm_external_ref_from_gsm_soft(gsm_soft, gsm):

    gsm_e_r = {}

    gsm_m = {
        "id": "!Sample_geo_accession = ",
        "label": "!Sample_title = "
    }

    for k, m in gsm_m.items():
        for gsm_l in gsm_soft:
            if m in gsm_l:
                gsm_l = re.sub(m, "", gsm_l)
                v = gsm_l.strip()
                gsm_e_r.update({k: v})
                continue

    if not "id" in gsm_e_r:
        print(gsm)
        print(gsm_soft)

    gsm_e_r.update({"reference":byc["filter_definitions"]["GEOseries"]["reference"]["root"]+gsm_e_r["id"]})
    gsm_e_r.update({"id":"geo:"+gsm_e_r["id"]})

    return gsm_e_r

################################################################################

def _gpl_external_ref_from_gse_soft(gpl_soft, gpl):

    gpl_e_r = {}

    gpl_m = {
        "id": "!Platform_geo_accession = ",
        "label": "!Platform_title = "
    }

    for k, m in gpl_m.items():
        for gpl_l in gpl_soft:
            if m in gpl_l:
                gpl_l = re.sub(m, "", gpl_l)
                v = gpl_l.strip()
                gpl_e_r.update({k: v})
                continue

    gpl_e_r.update({"reference":byc["filter_definitions"]["GEOseries"]["reference"]["root"]+gpl_e_r["id"]})
    gpl_e_r.update({"id":"geo:"+gpl_e_r["id"]})

    return gpl_e_r

################################################################################

def _read_retrieve_save_gse_soft(gse, byc):

    i_d = re.sub("geo:", "", gse)

    gse_path = path.join( dir_path, *byc["geosoft_file_root"], "GSE", i_d )
    gse_file_path = path.join( gse_path, i_d+".txt" )

    if path.isfile(gse_file_path):               
        gse_soft = open(gse_file_path).read().splitlines()
        return gse_soft, 0
    else:
        gse_soft = retrieve_geosoft_file(i_d, byc)
        if not path.isdir(gse_path):
            mkdir(gse_path)
        s_f = open(gse_file_path, 'w')
        for l in gse_soft:
            s_f.write(l)
        s_f.close()
        return gse_soft, 1

    return False, 0

################################################################################

def _read_retrieve_save_gpl_soft(gpl, byc):

    i_d = re.sub("geo:", "", gpl)

    gpl_path = path.join( dir_path, *byc["geosoft_file_root"], "GPL" )
    gpl_file_path = path.join( gpl_path, i_d+".txt" )

    if path.isfile(gpl_file_path):               
        gpl_soft = open(gpl_file_path).read().splitlines()
        return gpl_soft, 0
    else:
        gpl_soft = retrieve_geosoft_file(i_d, byc)
        if not path.isdir(gpl_path):
            mkdir(gpl_path)
        s_f = open(gpl_file_path, 'w')
        for l in gpl_soft:
            s_f.write(l)
        s_f.close()
        return gpl_soft, 1

    return False, 0

################################################################################

def _filter_gsm_soft(gsm_soft):

    # getting rid of data header and some verbose stuff
    gsm_soft = list(filter(lambda x:'Sample_' in x, gsm_soft))
    gsm_soft = list(filter(lambda x:'_protocol' not in x, gsm_soft))
    
    return gsm_soft

################################################################################

def _save_tmp_file(filename, content, byc):

    f_p = path.join( dir_path, *byc["config"]["paths"]["tmp_file_root"], filename )
    t_f = open(f_p, 'w')
    for l in content:
        t_f.write(l+"\n")
    t_f.close()

    return f_p

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
