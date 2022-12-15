import argparse, base36, csv, re, requests
from random import sample as randomSamples

from args_parsing import create_args_parser

################################################################################

def get_args(byc):

    set_processing_modes(byc)
    filters_from_args(byc)

    return byc

################################################################################

def set_single_dataset(byc):

    if len(byc["dataset_ids"]) > 1:
        print("Please give only one dataset using -d")
        exit()

    return byc["dataset_ids"][0]

################################################################################

def set_processing_modes(byc):

    byc.update({"update_mode": False})

    try:
        if byc["test_mode"] is True:
            byc.update({"update_mode": False})
            print( "¡¡¡ TEST MODE - no db update !!!")
            return byc
    except:
        pass

    try:
        if byc["args"].update:
            byc.update({"update_mode": True})
            print( "¡¡¡ UPDATE MODE - may overwrite entries !!!")
    except:
        pass
        
    return byc

################################################################################

def filters_from_args(byc):

    if not "args" in byc:
        return

    if not "filters" in byc:
        byc.update({"filters":[]})

    if byc["args"].filters:
        for f in re.split(",", byc["args"].filters):
            byc["filters"].append({"id":f})
 
    return byc

################################################################################

def genome_binning_from_args(byc):

    if byc["args"].key:
        byc.update({"genome_binning": byc["args"].key})

    return byc

