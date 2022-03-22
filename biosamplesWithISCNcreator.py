#!/usr/local/bin/python3

import re, json, yaml, statistics, sys, datetime, time
from os import path, environ, pardir
from pymongo import MongoClient
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

    biosamples_ISCN_updater()

################################################################################

def biosamples_ISCN_updater():

    initialize_service(byc, "biosamples_refresher")
    get_args(byc)
    set_processing_modes(byc)

    if not byc["args"].inputfile:
        print("No inputfile file specified => quitting ...")
        exit()

    ds_id = "progenetix"
    parse_filters(byc)
    parse_variants(byc)
    initialize_beacon_queries(byc)
    generate_genomic_intervals(byc)

    io_params = byc["datatable_mappings"]["io_params"]
    io_prefixes = byc["datatable_mappings"]["io_prefixes"]

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

    max_count = 1000000

    bios_temp = object_instance_from_schema_name(byc, "pgxBiosample")

    logfile = path.splitext(byc["args"].inputfile)[0]
    logfile += "_log.tsv"
    log = open(logfile, "w")

    fmp_samples, fmp_pmids, fieldnames = _read_samplefile(byc, max_count)
    fmp_no = len(fmp_samples)
    
    bar = Bar("{} samples from {} studies".format(fmp_no, len(fmp_pmids)), max = fmp_no, suffix='%(percent)d%%'+" of "+str(fmp_no) )

    for f_s in fmp_samples:

        bar.next()

        exists = False

        if "pgx" in f_s["biosample_id"]:
            exists = True
        else:
            b_id = legacy_ids.get(f_s["pgx_legacy_sample_id"], None)
            if b_id is not None:
                f_s.update({"biosample_id":b_id})
                print("found existing {} for {}".format(f_s["biosample_id"], f_s["pgx_legacy_sample_id"]))
                exists = True

        if exists is True:
            bios = bios_coll.find_one({"id": f_s["biosample_id"]})
            bios_id = bios["id"]
            ind_id = bios["individual_id"]
            cs = cs_coll.find_one({"biosample_id": bios_id})
            cs_id  = cs["id"]
            ind = ind_coll.find_one({"id": ind_id})
        else:
            bios_id = generate_id("pgxbs")
            bs_legacy_ids.update({bios_id: f_s["pgx_legacy_sample_id"] })
            f_s.update({"biosample_id":bios_id})
            cs_id = re.sub("pgxbs", "pgxcs", bios_id)
            ind_id = re.sub("pgxbs", "pgxind", bios_id)
            bios = deepcopy(bios_temp)
            bios.update({
                "id": bios_id,
                "individual_id": ind_id,
                "callset_ids": [ cs_id ],
                "provenance_note": "FMP recovery",
                "updated": isoformat(datetime.datetime.now())
            })
            bios["info"].update({
                "legacy_ids": [f_s["pgx_legacy_sample_id"]],
                "provenance_note": "FMP recovery"
            })
            ind = {
                "id": ind_id,
                "info":{"provenance_note": "FMP recovery"},
                "updated": isoformat(datetime.datetime.now())
            }
            cs = {
                "id": cs_id,
                "biosample_id": bios_id,
                "individual_id": ind_id,
                "info":{"provenance_note": "FMP recovery"},
                "updated": isoformat(datetime.datetime.now())
            }

        if not "biosample_id" in f_s:
            print(f_s)
            exit()

        bios = import_datatable_dict_line(byc, bios, fieldnames, f_s, "biosample")
        cs = import_datatable_dict_line(byc, cs, fieldnames, f_s, "analysis")
        ind = import_datatable_dict_line(byc, ind, fieldnames, f_s, "individual")

        _recreate_all_variants_from_ISCN(f_s, cs, v_coll, byc)

        if not byc["test_mode"]:
            bios_coll.update_one({"id": bios_id}, { "$set": bios },  upsert=True)
            ind_coll.update_one({"id": ind_id}, { "$set": ind },  upsert=True)
            cs_coll.update_one({"id": cs_id}, { "$set": cs },  upsert=True)

    bar.finish()

    if not byc["test_mode"]:
        with open(path.join(byconeer_path, *byc["this_config"]["biosample_id_path"]), 'w') as f:
            yaml.dump(bs_legacy_ids, f)

    exit()

################################################################################
##### TBD but probably elsewhere ###############################################
################################################################################

    # WARNING - ONLY FOR THE NEW SAMPLES

def _recreate_all_variants_from_ISCN(sample, cs, v_coll, byc):

    if not "ISCN" in sample["data_provenance"]:
        print("¡¡¡ No ISCN in data_provenance - skipping !!!")
        return

    cs_id = cs["id"]
    bs_id = cs["biosample_id"]

    technique = "aCGH"
    iscn_field = "iscn_acgh"
    platform_id = "EFO:0002701"
    platform_label = "DNA array"

    if "cCGH" in sample["data_provenance"]:
        technique = "cCGH"
        iscn_field = "iscn_ccgh"
        platform_id = "EFO:0010937"
        platform_label = "comparative genomic hybridization (CGH)"

    cs.update({ "platform_model": {"id": platform_id, "label":platform_label} })
    cs["info"].update( { "provenance": technique+" ISCN conversion" } )

    # print("\n{} ({}): {} - {}".format(bs_id, cs_id, technique, sample[iscn_field]))
    variants, variant_error = _deparse_rev_ish_CGH(bs_id, cs_id, technique, sample[iscn_field])

    if not byc["test_mode"]:

        v_coll.delete_many({"biosample_id": bs_id })

        for v in variants:
            v["info"].update({"provenance_note": "FMP recovery"})
            v_id = v_coll.insert_one(v).inserted_id
            v_coll.update_one({"_id": v_id}, {"$set": {"id": str(v_id)}})

        maps, cs_cnv_stats, cs_chro_stats = interval_cnv_arrays(v_coll, { "callset_id": cs_id }, byc)
        cs.update({"cnv_statusmaps": maps})
        cs.update({"cnv_stats": cs_cnv_stats})
        cs.update({"cnv_chro_stats": cs_chro_stats})

    return cs

################################################################################

def _read_samplefile(byc, max_count=0):

    fmp_samples = []
    fmp_pmids = set()
    fieldnames = []

    with open(byc["args"].inputfile, newline='') as csvfile:
        
        fmp_in = csv.DictReader(csvfile, delimiter="\t", quotechar='"')

        fieldnames += fmp_in.fieldnames

        for fmp_s in fmp_in:

            fmp_s = dict(fmp_s)

            pmid = fmp_s["external_references__id___PMID"]
            if not "PMID:" in pmid:
                pmid = "PMID:"+pmid

            fmp_s.update({"external_references__id___PMID":pmid})

            fmp_pmids.add(pmid)
            fmp_samples.append(dict(fmp_s))

    if max_count >0:
        if max_count < len(fmp_samples):
            fmp_samples = randomSamples(fmp_samples, k=max_count)

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

def _deparse_rev_ish_CGH(bs_id, cs_id, technique, iscn):

    v_s, v_e = deparse_ISCN_to_variants(iscn, technique, byc)
    variants = []

    for v in v_s:

        v.update({
            "variant_internal_id": variant_create_digest(v),
            "biosample_id": bs_id,
            "callset_id": cs_id,
            "updated": datetime.datetime.now().isoformat()
        })

        variants.append(v)

    return variants, v_e

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    main()
