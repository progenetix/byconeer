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
    _get_args(byc)

    byc["dataset_ids"] = ["progenetix"]

    parse_filters(byc)
    parse_variant_parameters(byc)
    initialize_beacon_queries(byc)
    generate_genomic_intervals(byc)

    print("Running progenetix ...")

    no_cs_no = 0
    no_stats_no = 0
    ds_id = byc["dataset_ids"][0]
    id_filter = "GSM"

    if byc["args"].filters:
        id_filter = byc["args"].filters

    bios_query = { "$or": [
        { "external_references.id": { "$regex": id_filter } },
        { "info.legacy_ids": { "$regex": id_filter } }
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

    missing_ids = []
    coll_lines = []
    coll_objs = []

    # header = ["id", "analysis_info.experiment_id", "analysis_info.series_id", "analysis_info.platform_id"]
    # for scp in sel_scopes:
    #     header.append(s_scopes[scp]["db_key"])
    #     header.append("_old_"+scp)
    #     header.append("_input_"+scp)
    #     header.append("_note_"+scp)

    # coll_lines.append("\t".join(header))

    for bsid in bs_ids:

        get_geosoft = True
        gsm = None
        gse = None
        gpl = None
        gsm_soft = False

        s = bios_coll.find_one({ "id":bsid })
        a_i = s.get("analysis_info", {})

        e_r_s = s.get("external_references", [])
        for e_r in e_r_s:
            if "GSM" in e_r["id"]:
                gsm = experiment_re.match(e_r["id"]).group(1)
            if "GSE" in e_r["id"]:
                print(e_r["id"])
                gse = series_re.match(e_r["id"]).group(1)

        if gsm is None:
            gsm = s["analysis_info"].get("experiment_id", None)
        if gse is None:
            gse = s["analysis_info"].get("series_id", None)

        if gsm is None:
            missing_ids.append("\t".join([bsid, gse, gsm, gpl]))
            bar.next()
            continue

        # changing this now to get GSE and GPL from GSM soft file

            # if "GSE" in e_r["id"]:
            #     gse = series_re.match(e_r["id"]).group(1)
            # if "GPL" in e_r["id"]:
            #     gpl = platform_re.match(e_r["id"]).group(1)

        gsm_soft = _read_retrieve_save_gsm_soft(gsm, gse, byc)

        if not gsm_soft:
            print("\n!!!! No soft file for {}".format(gsm))
            bar.next()
            continue

        if gse is None:
            gse = _return_gse_from_gsm_soft_list(gsm_soft)
        if gpl is None:
            gpl = _return_gpl_from_gsm_soft_list(gsm_soft)
        gse_soft = _read_retrieve_save_gse_soft(gse, byc)

        _gse_external_ref_from_gse_soft(gse_soft)

        _filter_gsm_soft(gsm_soft)

        updating_scopes = False

        line_coll = {}
        # for h_i in header:
        #     line_coll.update({h_i:""})
        # line_coll.update({
        #     "id": bsid,
        #     "analysis_info.experiment_id": gsm,
        #     "analysis_info.series_id": gse,
        #     "analysis_info.platform_id": gpl
        # })

        new_info = s["info"].copy()

        # sample_characteristics = geosoft_preclean_sample_characteristics(gsm_soft, byc)

        # coll_line = ""
        # if updating_scopes:
        #     line = []
        #     for h_k in header:
        #         line.append(line_coll.get(h_k, ""))
        #     coll_line = "\t".join(line)
        #     coll_lines.append(coll_line)
        #     coll_objs.append(line_coll)
        #     bar.next()
        # else:
        #     bar.next()
        #     continue

    bar.finish()

    # tmp_path = _save_tmp_file("gsm-metadata_"+"_".join(sel_scopes)+".tsv", coll_lines, byc)
    # print("=> Wrote {}".format(tmp_path))
    # print("=> metadata for {} samples".format(len(coll_lines) - 1))

    # for scp in sel_scopes:
    #     scp_dists = {}
    #     for c in coll_objs:
    #         if len(c[ s_scopes[scp]["db_key"] ]) > 0:
    #             scp_dists.update({ c[ s_scopes[scp]["db_key"] ] :1})
    #     print("=> Values in scope \"{}\":\n{}".format(s_scopes[scp]["id"], "\n".join(list(scp_dists.keys()))))

################################################################################

def _read_retrieve_save_gsm_soft(gsm, gse, byc):

    gsm_soft = None

    if gse is None:
        gsm_soft = retrieve_geosoft_file(gsm, byc)
        gse = _return_gse_from_gsm_soft_list(gsm_soft)

    gse_path = path.join( dir_path, *byc["geosoft_file_root"], gse )
    gse_file_path = path.join( gse_path, gse+".txt" )
    gsm_file_path = path.join( gse_path, gsm+".txt" )

    if path.isfile(gsm_file_path):
        if gsm_soft is None:            
            gsm_soft = open(gsm_file_path).read().splitlines()
        return gsm_soft
    else:
        if gsm_soft is None:
            gsm_soft = retrieve_geosoft_file(gsm, byc)
        if not path.isdir(gse_path):
            mkdir(gse_path)
        s_f = open(gsm_file_path, 'w')
        for l in gsm_soft:
            s_f.write(l)
        s_f.close()
        return gsm_soft

    return False

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

def _gse_external_ref_from_gse_soft(gse_soft):

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

    gse_e_r.update({"id":"geo:"+gse_e_r["id"]})

    prjsoncam(gse_e_r)
    return gse_e_r


################################################################################

def _read_retrieve_save_gse_soft(gse, byc):

    gse_path = path.join( dir_path, *byc["geosoft_file_root"], gse )
    gse_file_path = path.join( gse_path, gse+".txt" )

    if path.isfile(gse_file_path):               
        gse_soft = open(gse_file_path).read().splitlines()
        return gse_soft
    else:
        gse_soft = retrieve_geosoft_file(gse, byc)
        if not path.isdir(gse_path):
            mkdir(gse_path)
        s_f = open(gse_file_path, 'w')
        for l in gse_soft:
            s_f.write(l)
        s_f.close()
        return gse_soft

    return False

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
