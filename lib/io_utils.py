import argparse

################################################################################

def get_args(byc):

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datasetids", help="datasets, comma-separated")
    parser.add_argument("-a", "--alldatasets", help="process all datasets", action='store_true' )
    parser.add_argument("-c", "--collationtypes", help='selected collation types, e.g. "EFO"')
    parser.add_argument("-f", "--filters", help="prefixed filter values, comma concatenated")
    parser.add_argument("-t", "--test", help="test setting")
    parser.add_argument('-i', '--inputfile', help='a custom file to specify input data')
    parser.add_argument("-o", "--outfile", help="output file")
    parser.add_argument('-m', '--mode', help='update modus')
    parser.add_argument("-s", "--source", help="id source", nargs='?', default="callsets")
    parser.add_argument("-u", "--update", help='update existing records')
    parser.add_argument("-y", "--ontologycodes", help='selected codes , e.g. "NCIT:C2955"')

    byc.update({ "args": parser.parse_args() })

    return byc

################################################################################

def set_test_mode(byc):

    byc.update({"test_mode": False})

    try:
        if byc["args"].test:
            byc.update({"test_mode": True})
            print( "¡¡¡ TEST MODE - no db update !!!")
    except:
        pass
        
    return byc

################################################################################
