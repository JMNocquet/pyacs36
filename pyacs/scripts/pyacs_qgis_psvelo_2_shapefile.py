#!/usr/bin/env python

###################################################################
# SCRIPT    : pyacs_qgis_psvelo_2_shapefile.py
# AUTHOR    : JM NOCQUET
# DATE      : September 2012
# INPUT     : gmt file
# OUTPUT    : shapefile
# NOTE      :
#         
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import os,sys
import argparse
import pyacs.lib.shapefile

def file_base_name(file_name):
    if '.' in file_name:
        separator_index = file_name.index('.')
        base_name = file_name[:separator_index]
        return base_name
    else:
        return file_name

def path_base_name(path):
    file_name = os.path.basename(path)
    return file_base_name(file_name)


###################################################################
# PARSE ARGUMENT LINE
###################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-gmt', action='store', type=str, dest='gmt',required=True,help='gmt psvelo file (long. lat. Ve Vn S_Ve S_Vn S_Ven name)')
parser.add_argument('-shapefile', action='store', type=str, dest='shapefile',default=None,help='output shapefile name')
parser.add_argument('--verbose', action='count', default=0, help='verbose mode')

if (len(sys.argv)<2):parser.print_help()
args = parser.parse_args()

if args.shapefile is None:
    shp = path_base_name( args.gmt )
else:
    shp = args.shapefile

if args.verbose>0:
    verbose=True
else:
    verbose=False

pyacs.lib.shapefile.psvelo_to_shapefile(args.gmt, shp, verbose=verbose)
