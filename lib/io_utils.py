import argparse, base36, csv, re, time
from random import sample as randomSamples

################################################################################

def get_args(byc):

    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alldatasets", help="process all datasets", action='store_true' )
    parser.add_argument("-c", "--collationtypes", help='selected collation types, e.g. "EFO"')
    parser.add_argument("-d", "--datasetids", help="datasets, comma-separated")
    parser.add_argument("-f", "--filters", help="prefixed filter values, comma concatenated")
    parser.add_argument('-i', '--inputfile', help='a custom file to specify input data')
    parser.add_argument('-m', '--mode', help='update modus')
    parser.add_argument("-o", "--outfile", help="output file")
    parser.add_argument("-q", "--query", help='complete query string, e.g. `{"biosamples":{"external_references.id":"geo:GSE7428"}}`')
    parser.add_argument("-r", "--randno", help="random number", default=0)
    parser.add_argument("-s", "--source", help="some source", nargs='?', default="callsets")
    parser.add_argument("-t", "--test", help="test setting")
    parser.add_argument("-u", "--update", help='update existing records')
    parser.add_argument("-y", "--ontologycodes", help='selected codes , e.g. "NCIT:C2955"')

    byc.update({ "args": parser.parse_args() })

    return byc

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

def set_processing_modes(byc):

    byc.update({
        "test_mode": False,
        "update_mode": False
    })

    try:
        if byc["args"].test:
            byc.update({"test_mode": True})
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

def generate_id(prefix):

    time.sleep(.001)
    return '{}-{}'.format(prefix, base36.dumps(int(time.time() * 1000))) ## for time in ms
