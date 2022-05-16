import argparse, base36, csv, re, time
from random import sample as randomSamples

from args_parsing import create_args_parser

################################################################################

def get_args(byc):

    # print(byc["pkg_path"])
    # parser.add_argument("-c", "--collationtypes", help='selected collation types, e.g. "EFO"')
    # parser.add_argument("-q", "--query", help='complete query string, e.g. `{"biosamples":{"external_references.id":"geo:GSE7428"}}`')

    # byc.update({ "args": parser.parse_args() })

    # if "script_args" in byc:
    #     create_args_parser(byc, byc["script_args"])

    set_processing_modes(byc)

    return byc

################################################################################

def set_single_dataset(byc):

    if len(byc["dataset_ids"]) > 1:
        print("Please give only one dataset using -d")
        exit()
    return byc["dataset_ids"][0]

################################################################################

def read_tsv_to_dictlist(filepath, max_count=0):

    dictlist = []

    with open(filepath, newline='') as csvfile:
    
        data = csv.DictReader(csvfile, delimiter="\t", quotechar='"')
        fieldnames = list(data.fieldnames)

        for l in data:
            dictlist.append(dict(l))

    if max_count > 0:
        if max_count < len(dictlist):
            dictlist = randomSamples(dictlist, k=max_count)

    return dictlist, fieldnames

################################################################################

def read_inputtable(byc, max_count=0):

    ind_data = []
    bs_ids = set()
    fieldnames = []

    with open(byc["args"].inputfile, newline='') as csvfile:
        
        fmp_in = csv.DictReader(csvfile, delimiter="\t", quotechar='"')

        fieldnames += fmp_in.fieldnames

        for fmp_s in fmp_in:

            fmp_s = dict(fmp_s)

            bsid = fmp_s["biosample_id"]
            bs_ids.add(bsid)
            ind_data.append(dict(fmp_s))

    if max_count >0:
        if max_count < len(ind_data):
            ind_data = randomSamples(ind_data, k=max_count)

    return ind_data, list(bs_ids), fieldnames

################################################################################

def set_processing_modes(byc):

    byc.update({"update_mode": False})

    try:
        if byc["test_mode"] is True:
            print( "¡¡¡ TEST MODE - no db update !!!")
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

    return(byc)

################################################################################

def generate_id(prefix):

    time.sleep(.001)
    return '{}-{}'.format(prefix, base36.dumps(int(time.time() * 1000))) ## for time in ms
