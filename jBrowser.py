#!/usr/local/bin/python3

import re, json, yaml
from os import path, environ, pardir
import sys, datetime
from pymongo import MongoClient

from dash import html, Dash
from jbrowse_jupyter import create, create_component

# bycon is supposed to be in the same parent directory
dir_path = path.dirname( path.abspath(__file__) )
parent_path = path.join( dir_path, pardir )
sys.path.append( parent_path )

from bycon import *
from byconeer import *

app = Dash(__name__)

jbrowse_conf = create("LGV", genome="hg38")

config = jbrowse_conf.get_config()

component = create_component(config)

app.layout = html.Div(
    [component],
    id='test'
)

################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
    app.run_server(port=8081, debug=True)
