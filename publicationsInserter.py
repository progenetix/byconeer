#!/usr/local/bin/python3

from os import path, pardir
from pymongo import MongoClient
from isodate import date_isoformat
import cgi, cgitb, csv, datetime, requests, sys

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""
* pubUpdater.py -t 1 -f "../rsrc/publications.txt"
"""

##############################################################################
##############################################################################
##############################################################################

def main():
    publications_inserter()

##############################################################################

def publications_inserter():

    initialize_service(byc)
    get_args(byc)
    set_processing_modes(byc)

    g_url = byc["service_config"]["google_spreadsheet_tsv_url"]
    skip_cols = byc["service_config"]["skipped_columns"]

    if byc["args"].inputfile:
        pub_file = yc["args"].inputfile
    else:
        print("No inputfile file specified => pulling the online table ...")
        # exit()
        pub_file = path.join( dir_path, "tmp", "pubtable.tsv" )
        print("... reading from {}".format(g_url["base_url"]))
        r =  requests.get(g_url["base_url"], params=g_url["params"])
        if r.ok:
            with open(pub_file, 'wb') as f:
                f.write(r.content)
            print("Wrote file to {}".format(pub_file))
        else:
            print("Download failed: status code {}\n{}".format(r.status_code, r.text))

    rows = []

    mongo_client = MongoClient()

    pub_coll = mongo_client["progenetix"]["publications"]
    bios_coll = mongo_client["progenetix"]["biosamples"]

    publication_ids = pub_coll.distinct("id")

    progenetix_ids = bios_coll.distinct("external_references.id")
    progenetix_ids = [item for item in progenetix_ids if item is not None]
    progenetix_ids = list(filter(lambda x: x.startswith("PMID"), progenetix_ids))

    # TODO: Use schema ...

    up_count = 0

    with open(pub_file, newline='') as csvfile:

        in_pubs = list(csv.DictReader(csvfile, delimiter="\t", quotechar='"'))

        print("=> {} publications will be looked up".format(len(in_pubs)))

        l_i = 0
        for pub in in_pubs:

            l_i += 1

            pmid = str(pub.get("pubmedid", "empty")).strip()
            skip_mark = pub.get("SKIP", "").strip()

            if len(skip_mark) > 0:
                print('¡¡¡ Line {} ({}): skipped due to non-empty skip field ("{}") !!!'.format(l_i, pmid, skip_mark))
                continue

            if not re.match(r'^\d{6,9}$', pmid):
                print('¡¡¡ Line {}: skipped due to empty or strange pubmedid entry ("{}") !!!'.format(l_i, pmid))
                continue

            p_k = "PMID:"+pmid

            """Publications are either created from an empty dummy or - if id exists and
            `-u 1` taken from the existing one."""

            if p_k in publication_ids:
                if not byc["update_mode"]:
                    print(p_k, ": skipped - already in progenetix.publications")
                    continue
                else:
                    n_p = mongo_client["progenetix"]["publications"].find_one({"id": p_k })
                    print(p_k, ": existed but overwritten since *update* in effect")
            else:
                n_p = get_empty_publication(byc)
                n_p.update({"id":p_k})

            for k, v in pub.items():
                v = v.strip()
                if k:
                    if k in skip_cols:
                        continue
                    if len(str(v)) < 1:
                        continue
                    if v == "DELETE":
                        v = ""
                    assign_nested_value(n_p, k, v)

            try:
                if len(pub["PROVENANCE_ID"]) > 4:
                    geo_info = mongo_client["progenetix"]["geolocs"].find_one({"id": pub["PROVENANCE_ID"]}, {"_id": 0, "id": 0})
                    if geo_info is not None:
                        n_p["provenance"].update({"geo_location":geo_info["geo_location"]})
            except KeyError:
                pass

            epmc, e = retrieve_epmc_publications(pmid)
            if e is not False:
                print(e)
                continue

            update_from_epmc_publication(n_p, epmc)            
            publication_update_label(n_p)
            get_ncit_tumor_types(n_p, pub)

            if p_k in progenetix_ids:

                n_p["counts"].update({ "progenetix" : 0 })
                n_p["counts"].update({ "arraymap" : 0 })

                for s in bios_coll.find({ "external_references.id" : p_k }):
                    n_p["counts"]["progenetix"] += 1
                for s in bios_coll.find({ "cohorts.id" : "pgxcohort-arraymap" }):
                    n_p["counts"]["arraymap"] += 1

            for c in n_p["counts"].keys():
                if isinstance(n_p["counts"][c], str):
                    try:
                        n_p["counts"].update({c: int(n_p["counts"][c])})
                    except:
                        pass
            n_p["counts"]["ngs"] = n_p["counts"]["wes"] + n_p["counts"]["wgs"]

            if not byc["test_mode"]:
                entry = pub_coll.update_one({"id": n_p["id"] }, {"$set": n_p }, upsert=True )
                up_count += 1
                print(n_p["id"]+": inserting this into progenetix.publications")
            else:
                jprint(n_p)
                    
    print("{} publications were inserted or updated".format(up_count))

##############################################################################
##############################################################################



##############################################################################
##############################################################################

if __name__ == '__main__':
        main()

##############################################################################
##############################################################################
##############################################################################
