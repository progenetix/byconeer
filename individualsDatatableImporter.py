#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir, mkdir
import sys, datetime
from pymongo import MongoClient
import statistics
from progress.bar import Bar
import requests

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

################################################################################
################################################################################
################################################################################

def main():

    individuals_refresher()

################################################################################

def individuals_refresher():

    initialize_service(byc)
    set_processing_modes(byc)
    select_dataset_ids(byc)

    if not byc["args"].inputfile:
        print("No inputfile file specified => quitting ...")
        exit()
    
    ds_id = byc["dataset_ids"][0]
    max_count = 1000000


    # ind_temp = object_instance_from_schema_name(byc, "individual", "properties")
    # print(ind_temp)

    data_client = MongoClient( )
    data_db = data_client[ ds_id ]
    bios_coll = data_db[ "biosamples" ]
    ind_coll = data_db[ "individuals" ]

    ind_table, bs_ids, fieldnames = read_inputtable(byc, max_count)
    ind_no = len(ind_table)

    for t_ind in ind_table:

        bs_id = t_ind.get("biosample_id", None).strip()

        if not bs_id:
            continue
        bs = bios_coll.find_one({"id": bs_id})
        if not bs:
            print("¡¡¡ biosample {} could not be found !!!".format(bs_id))
            continue
        ind_id = bs["individual_id"]
        ind = ind_coll.find_one({"id":ind_id})
        if not ind:
            print("¡¡¡ individual {} could not be found !!!".format(ind_id))
            continue

        update_bs = import_datatable_dict_line(byc, {}, fieldnames, t_ind, "biosample")
        update_ind = import_datatable_dict_line(byc, {}, fieldnames, t_ind, "individual")
        if byc["update_mode"] is True:

            if len(update_ind.keys()) > 0:
                ind_coll.update_one({"_id": ind["_id"] }, { "$set": update_ind },  upsert=False)

            if len(update_bs.keys()) > 0:
                bios_coll.update_one({"_id": bs["_id"] }, { "$set": update_bs },  upsert=False)

        else:
            prjsonnice(update_ind)
            prjsonnice(update_bs)

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
