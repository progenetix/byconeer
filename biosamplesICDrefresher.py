#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from pymongo import MongoClient
import statistics
from progress.bar import Bar

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""

## `biosamplesRefresher`

* `biosamplesICDrefresher.py -d progenetix

"""

################################################################################
################################################################################
################################################################################

ds_id = "progenetix"
data_client = MongoClient( )
data_db = data_client[ ds_id ]
bios_coll = data_db[ "biosamples" ]

pgx_sample_link = "https://progenetix.org/services/ids/"

def main():

    biosamples_refresher()

################################################################################

def biosamples_refresher():

    initialize_service(byc)
    parse_filters(byc)
    parse_variant_parameters(byc)
    initialize_beacon_queries(byc)

    icdo_ncit_mapping_file = path.join( dir_path, "rsrc", "NCIT", "icdom_ncit_progenetix_mappings.tsv")

    modes = {
        "icdom": {
            "input_file": path.join( dir_path, "rsrc", "icdom", "ICD-O-3.1-NCIt_Morphology_Mapping.tsv"),
            "parameter": "icdo_morphology",
            "pgx_hierarchy_file": path.join( dir_path, "rsrc", "icdom", "numbered_hierarchies.tsv"),
            "update_no": 0
        },
        "icdot": {
            "input_file": path.join( dir_path, "rsrc", "icdot", "ICD-O-3.1-NCIt_Topography_Mapping.tsv"),
            "parameter": "icdo_topography",
            "pgx_hierarchy_file": path.join( dir_path, "rsrc", "icdot", "numbered_hierarchies.tsv"),
            "update_no": 0
        },
    }

    for m_k, m_p in modes.items():
        
        m_p.update({"codes": _map_codes_from_file(m_p["input_file"], m_k) })
        id_p = m_p["parameter"]+".id"
        label_p = m_p["parameter"]+".label"

        no = len(m_p["codes"])

        if not byc["test_mode"]:
            bar = Bar("{} {} codes".format(no, m_k), max = no, suffix='%(percent)d%%'+" of "+str(no) )

        for icd_k, icd_l in m_p["codes"].items():

            q = {id_p: icd_k}

            f_l = bios_coll.count_documents(q)


            if f_l > 0:
                m_p.update({ "update_no": m_p["update_no"] + f_l })
                if not byc["test_mode"]:
                    bios_coll.update_many(q, {"$set":{label_p:icd_l}})
                else:
                    print("{}: {} ({})".format(icd_k, icd_l, f_l))

            if not byc["test_mode"]:
                bar.next()

        if not byc["test_mode"]:
            bar.finish()

        print("{} samples updated for {}".format(m_p["update_no"], label_p))

        _rewrite_hierarchy_file(m_p["input_file"], m_p["pgx_hierarchy_file"], m_k)
        
    _rewrite_icdo_ncit_mapping_file(path.join( dir_path, "rsrc", "icdom", "ICD-O-3.1-NCIt_Morphology_Mapping.tsv"), icdo_ncit_mapping_file, modes)

################################################################################

def _map_codes_from_file(input_file, mode):

    mapped = {}

    with open(input_file) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('#'):
                continue
            # ICD-O Code    Level   Term Type   Term Type Desc  ICD-O string    NCIm CUI    NCIt Code (if present)  NCIt PT string (Preferred term)
            icd_c, level, t_type, t_desc, icd_l, ncit_cui, ncit_c, ncit_l = line.split("\t")
            if not "PT" in t_type:
                continue

            if "icdom" in mode:
                icd_k = "pgx:icdom-{}".format(re.sub('/', '', icd_c))
            elif "icdot" in mode:
                icd_k = "pgx:icdot-{}".format(icd_c)

            mapped[icd_k] = icd_l

    return mapped

################################################################################

def _rewrite_hierarchy_file(input_file, pgx_hierarchy_file, mode):

    h_lines = []
    l_no = 0

    with open(input_file) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('#'):
                continue
            # ICD-O Code    Level   Term Type   Term Type Desc  ICD-O string    NCIm CUI    NCIt Code (if present)  NCIt PT string (Preferred term)
            icd_c, level, t_type, t_desc, icd_l, ncit_cui, ncit_c, ncit_l = line.split("\t")
            if "SY" in t_type:
                continue

            if not re.match(r'^\d$', level):
                continue

            level = int(level)

            if level == 1:
                continue

            level -= 2

            l_no += 1

            if "icdom" in mode:
                icd_k = "pgx:icdom-{}".format(re.sub('/', '', icd_c))
            elif "icdot" in mode:
                icd_k = "pgx:icdot-{}".format(icd_c)

            h_lines.append("{}\t{}\t{}\t{}".format(icd_k, icd_l, level, l_no))

    h_f = open(pgx_hierarchy_file, 'w')
    for l in h_lines:
        h_f.write(l+"\n")
    h_f.close()

################################################################################

def _rewrite_icdo_ncit_mapping_file(input_file, icdo_ncit_mapping_file, modes):

    h_lines = []
    l_no = 0

    mode = "icdom"

    id_p = modes[mode]["parameter"]+".id"
    icdoms = modes[mode]["codes"]

    h_lines.append("mapping_type\ticd_code_mix\ticdom_id\ticdom_label\ticdot_id\ticdot_label\tncit_id\tncit_label\tpgx_sample_example\tpgx_example_link\tpgx_sample_count")

    with open(input_file) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('#'):
                continue
            # ICD-O Code    Level   Term Type   Term Type Desc  ICD-O string    NCIm CUI    NCIt Code (if present)  NCIt PT string (Preferred term)
            icd_c, level, t_type, t_desc, icd_l, ncit_cui, ncit_c, ncit_l = line.strip("\n").split("\t")
            if not "PT" in t_type:
                continue

            level = 1
            l_no += 1
            icd_k = "pgx:icdom-{}".format(re.sub('/', '', icd_c))
            q = {id_p: icd_k}
            i_no = bios_coll.count_documents(q)

            if "C" in ncit_c:
                ncit_c = "NCIT:"+ncit_c
            h_lines.append("primary ICD-O M\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(re.sub("pgx:", "",icd_k), icd_k, icd_l, "", "", ncit_c, ncit_l, "", "", i_no))

            ncits = {}

            for s in bios_coll.find(q):
                k = "{}::{}".format(s["icdo_topography"]["id"], s["histological_diagnosis"]["id"])
                ncits[k] = {
                    "icdot": s["icdo_topography"],
                    "ncit": s["histological_diagnosis"],
                    "example_description": s.get("description", ""),
                    "example_id": pgx_sample_link + s.get("id", "")
                }

            for i_k in sorted(ncits.keys()):
                i_v = ncits[i_k]
                c_q = { id_p: icd_k, "icdo_topography.id":i_v["icdot"]["id"], "histological_diagnosis.id":i_v["ncit"]["id"] }
                c_no = bios_coll.count_documents(c_q)
                h_lines.append("pgx code mix\t{}~{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(re.sub("pgx:", "", icd_k), re.sub("pgx:", "", i_v["icdot"]["id"]), icd_k, icd_l, i_v["icdot"]["id"], i_v["icdot"]["label"], i_v["ncit"]["id"], i_v["ncit"]["label"], i_v["example_description"], i_v["example_id"], c_no))

    h_f = open(icdo_ncit_mapping_file, 'w')
    for l in h_lines:
        h_f.write(l+"\n")
    h_f.close()

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
