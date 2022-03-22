#!/usr/local/bin/python3
from pymongo import MongoClient
from os import path, environ, pardir
import sys, datetime, time
from progress.bar import Bar
from copy import deepcopy as deepcopy

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

"""
"""

################################################################################

def main():

	schema_parser()

################################################################################

def schema_parser():

	initialize_service(byc, "callsets_modifier")
	get_args(byc)
	if not byc["args"].parse:
		print("No schema name was provided trough the `-p` argument ...")
		exit()

	root_key = ""
	if byc["args"].key:
		root_key = byc["args"].key

	schema_name = byc["args"].parse

	s_f_p = get_schema_file_path(schema_name, byc)
	print("=> Using schema {} from {}".format(schema_name, s_f_p))
	s = read_schema_file(schema_name, root_key, byc)
	s_i = object_instance_from_schema_name(byc, schema_name)

	prjsonnice(s)
	print("\n################################################\n")
	prjsonnice(s_i)

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':
	main()
