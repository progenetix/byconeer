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

* `biosamplesICDrefresher.py

"""

################################################################################
################################################################################
################################################################################

def main():

    templates_creator()

################################################################################

def templates_creator():

    initialize_service(byc)
    dt_m = byc["datatable_mappings"]
    b_rt_s = byc["beacon_mappings"]["response_types"]

    rsrc_p = path.join(parent_path, "byconaut", "rsrc", "templates")

    added = []

    if byc["args"].parse:
        added = re.split(",", byc["args"].parse)

    for r_t in ["biosample", "individual", "variant", "analysis"]:

        io_params = dt_m["io_params"][ r_t ]
        io_prefixes = dt_m["io_prefixes"][ r_t ]

        coll = b_rt_s["biosample"].get("collection")

        header = create_table_header(io_params, io_prefixes)
        if len(added) > 0:
            header[1:1] = added

        f_p = path.join(rsrc_p, r_t+"s_template.tsv")
        f_p = re.sub("analysiss", "analyses", f_p)
        f = open(f_p, "w")
        f.write("\t".join(header)+"\n")
        f.close()

        k_s = list(io_params.keys())
        for io_p in k_s:
            d = io_params[io_p].get("compact", False)
            if d is False:
                io_params.pop(io_p, None)
        k_s = list(io_prefixes.keys())
        for io_p in k_s:
            d = io_prefixes[io_p].get("compact", False)
            if d is False:
                io_prefixes.pop(io_p, None)
            else:
                p_s = list(io_prefixes[io_p]["pres"].keys())
                for io_pre in p_s:
                    d_pre = io_prefixes[io_p]["pres"][io_pre].get("compact", False)
                    if d_pre is False:
                        io_prefixes[io_p]["pres"].pop(io_pre, None)

        header = create_table_header(io_params, io_prefixes)
        if len(added) > 0:
            header[1:1] = added

        f_p = path.join(rsrc_p, coll+"_compact_template.tsv")
        f = open(f_p, "w")
        f.write("\t".join(header)+"\n")
        f.close()

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
    main()
