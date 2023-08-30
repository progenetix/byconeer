#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from pymongo import MongoClient
from progress.bar import Bar
from pyexcel import get_sheet

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
pkg_path = path.join( dir_path, pardir )
parent_path = path.join( pkg_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""

## `ontologymapsCreator`

"""

################################################################################
################################################################################
################################################################################

def main():

    ontologymaps_creator()

################################################################################

def ontologymaps_creator():

    initialize_bycon_service(byc)
    run_beacon_init_stack(byc)

    if len(byc["dataset_ids"]) > 1:
        print("Please give only one dataset using -d")
        exit()

    ds_id = byc["dataset_ids"][0]
    print("=> Using data values from {}".format(ds_id))

    execute_bycon_queries( ds_id, byc )
    ds_results = byc["dataset_results"][ds_id]

    mongo_client = MongoClient( )
    data_db = mongo_client[ ds_id ]
    bios_coll = data_db[ "biosamples" ]

    if not "biosamples._id" in ds_results.keys():
        bios_ids = bios_coll.distinct("_id", {})
        print("¡¡¡ Using all {} biosamples from {} !!!".format(len(bios_ids), ds_id))
    else:
        bios_ids = ds_results["biosamples._id"]["target_values"]

    no = len(bios_ids)

    o_m_all = { }

    for mt, mv in byc["service_config"]["map_types"].items():

        o_m = { }

        o_l_max = len(mv["ontology_types"])
        
        if not byc["test_mode"]:
            bar = Bar("Analyzing samples", max = no, suffix="%(percent)d%%"+" of "+str(no) )

        for bios__id in bios_ids:

            s = bios_coll.find_one({"_id": bios__id})
            o_l_c = 0
            k_l = [ ]
            o_l = [ ]
            d = ""
            if "description" in s:
                d = s["description"].strip()

            o_t_l = len(mv["ontology_types"])

            for o_t in mv["ontology_types"]:

                data_key = byc["filter_definitions"][ o_t ]["db_key"]
                parent_key = re.sub(".id", "", data_key)
                data_re = re.compile( byc["filter_definitions"][ o_t ]["pattern"] )

                try:
                    o_p = s[ parent_key ]
                except:
                    continue

                if not o_p:
                    continue

                if not "id" in o_p:
                    continue

                if not data_re.match( o_p["id"] ):
                    continue

                if not "label" in o_p:
                    o_p.update({"label":""})

                k_l.append( o_p["id"] )
                o_l.append( o_p )

            if len(k_l) < o_t_l:
                continue
            k = "::".join(k_l)

            if k in o_m:
                if len(o_m[k]["examples"]) < byc["service_config"]["example_length"]:
                    o_m[k]["examples"].add(d)
                    o_m_all[k]["examples"].add(d)

            else:
                if byc["test_mode"]:
                    print(k)
                o_m[k] = {
                    "id": k,
                    "examples": set([d]),
                    "code_group": o_l,
                    "note": ds_id
                }

                o_m_all[k] = o_m[k].copy()
                o_m_all[k].pop("note")

            if byc["args"].filters:
                o_m[k].update({"note": byc["args"].filters})
                if "TCGA" in byc["args"].filters:
                    e_r_s = s.get("external_references", [])
                    for e_r in e_r_s:
                        if "pgx:TCGA" in e_r["id"]:
                            if "project" in e_r["label"]:
                                n = re.sub("pgx:", "", e_r["id"])
                                o_m[k].update({"note":n})

            if not byc["test_mode"]:
                bar.next()
        if not byc["test_mode"]:
            bar.finish()

        mt_l = mt

        if byc["args"].filters:
            mt_l = "{}-{}".format(mt, byc["args"].filters)
            mt_l = re.sub(":", "", mt_l)

        print("{} code combinations for {}".format(len(o_m.keys()), mt))
        _export_currentmaps(dir_path, mt_l, o_m)

        def_m = _read_mapping_defaults(dir_path, mt, **byc)

        for def_k in def_m:
            if not def_k in o_m.keys():
                o_m.update( { def_k: def_m[def_k] } )
        print("Now {} code combinations after defaults ...".format(len(o_m.keys())))
        _export_ontologymaps(dir_path, mt_l, o_m)

    if not byc["test_mode"]:

        for om_k, om_v in o_m_all.items():
            o_m_all[om_k].update({"examples": list(om_v["examples"]) })


        if not byc["args"].filters:
            om_coll = mongo_client["progenetix"]["ontologymaps"]
            om_coll.drop()
            om_coll.insert_many( o_m_all.values() )
            print("==> Rewrote {}.{} collection".format(byc["config"]["services_db"], byc["config"]["ontologymaps_coll"]))

################################################################################
################################################################################
################################################################################

def _export_currentmaps(dir_path, map_type, o_m):

    f_p = path.join(dir_path, "exports", "ontologymaps", "{}.tsv".format(map_type))
    f = open(f_p, "w")
    for o_k in sorted(o_m.keys()):
        l = [o_m[o_k]["note"]]
        for g_m in o_m[o_k]["code_group"]:
            l.append(str(g_m["id"]))
            l.append(str(g_m["label"]))
        l.append(str(sorted(o_m[o_k]["examples"])[0]))
        f.write("\t".join(l)+"\n")
    f.close()

################################################################################

def _export_ontologymaps(dir_path, map_type, o_m):
    for k, content in o_m.items():
        if "examples" in content:
            content["examples"] = sorted(content["examples"])
    export = [ ]
    for o_k in sorted(o_m.keys()):
        export.append(o_m[o_k])
    yaml.dump(export, open(path.join(dir_path, "exports", "ontologymaps", "{}.yaml".format(map_type)),"w"))

################################################################################
################################################################################
################################################################################

def _read_mapping_defaults(dir_path, map_type, **byc):

    if not "mappingfile" in byc["service_config"]["map_types"][ map_type ]:
        return {}

    mf = path.join( dir_path, *byc["service_config"]["map_types"][ map_type ]["mappingfile"] )
    o_m_r = { }
    equiv_keys = [ ]
    pre_fs = byc["service_config"]["map_types"][ map_type ]["ontology_types"]
    o_l_max = len(pre_fs)

    for o_t in pre_fs:
        equiv_keys.append( o_t+"::id" )
        equiv_keys.append( o_t+"::label" )

    sheet_name = "__".join(pre_fs)+"__matched"

    try:
        table = get_sheet(file_name=mf, sheet_name=sheet_name)
    except Exception as e:
        print(e)
        print("No matching mapping file could be found!")
        exit()

    header = table[0]
    col_inds = { }
    hi = 0
    fi = 0
    for col_name in header:
        if col_name in equiv_keys:
            col_inds[ col_name ] = hi

        hi += 1

    for i in range(1, len(table)):
        id_s = [ ]
        bioc_s = [ ]
        bioc = { }
        col_match_count = 0
        for col_name in equiv_keys:
            try:
                cell_val = table[ i, col_inds[ col_name ] ]
                if "id" in col_name:
                    o_t, code = re.split("[:-]", cell_val)
                    data_re = re.compile( byc["filter_definitions"][ o_t ]["pattern"] )
                    if data_re.match( cell_val ):
                        bioc = { "id": cell_val }
                        id_s.append( cell_val )
                else:
                    bioc.update( { "label": cell_val } )
                    bioc_s.append(bioc)
                    if len(id_s) == o_l_max:
                        o_k = "::".join(id_s)
                        o_m_r.update(
                            { o_k:
                                {
                                    "id": o_k,
                                    "code_group": bioc_s
                                }
                            }
                        )
                        fi += 1
            except:
                continue

    print("default mappings: "+str(fi))
    return o_m_r

################################################################################
################################################################################
################################################################################

if __name__ == "__main__":
    main()
