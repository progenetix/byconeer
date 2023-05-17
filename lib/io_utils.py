import argparse, base36, csv, re, requests
from random import sample as randomSamples

from args_parsing import *

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


