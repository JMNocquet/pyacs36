#!/usr/bin/env python

###################################################################
# SCRIPT    : pyacs_gvel_pole.py
# AUTHOR    : JM NOCQUET
# DATE      : September 2011 / April 2019
# INPUT     : gmt psvelo file
#           : experiment name
#           | input plate file
# OUTPUT    : 
# NOTE      :
#         
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import sys
import argparse

import numpy as np
from pyacs.lib.vel_field import Velocity_Field as VF
from pyacs.lib.errors import PyacsError

###################################################################
# PARSE ARGUMENT LINE
###################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-euler_file', action='store', type=str, dest='euler_file', \
                    help='file including plate name and euler pole values e.g. soam_itrf2008_nocquet_2014:-132.21/-18.83/0.121')
parser.add_argument('-plate', action='store', type=str, dest='plate', help='plate name from euler_file')
parser.add_argument('-value', action='store', type=str, dest='value', help='euler pole /long/lat/w in dec.deg / Myr')
parser.add_argument('-ivel', action='store', type=str, dest='ivel', help="input velocity file")
parser.add_argument('-ovel', action='store', type=str, dest='ovel', help="output velocity file. print if not provided")
parser.add_argument('-pll', action='append', type=str, dest='pll', help="show pole velocity prediction at /long/lat")
parser.add_argument('--verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('--debug', action='count',default=0,help='debug mode')

args = parser.parse_args()

# verbose
if args.verbose > 0:
    print("-- Verbose mode")
    verbose = True
else:
    verbose = False

# debug
if args.debug > 0:
    print("-- Debug mode")
    verbose = True
    debug = True
else:
    debug = False

###################################################################
# CHECK ARGUMENTS COMPATIBILITY
###################################################################

if (args.euler_file is None) and (args.plate is not None):
    print("!!! ERROR: argument plate provided but not euler_file")
    sys.exit()
    
if (args.value is not None) and (args.plate is not None):
    print("!!! ERROR: argument value and plate provided. PYACS does not know which one to choose.")
    sys.exit()

if (args.ivel is not None) and (args.pll is not None):
    print("!!! ERROR: argument ivel and pll provided. PYACS does not know which one to choose.")
    sys.exit()

if (args.ivel is None) and (args.pll is None):
    print("!!! ERROR: argument ivel or pll missing.")
    sys.exit()

###################################################################
# READ EULER FILE
###################################################################
H_euler={}
if args.euler_file is not None:
    print("-- reading Euler file: %s" % args.euler_file)
    try:
        ef=open(args.euler_file)
        for line in ef:
            if line[0]=='#':continue
            
            lline = line.split(':')
            name = lline[0].strip()
            value = np.array( list(map(float,lline[1].split('/') ) ) )
            H_euler[name] = value
        ef.close()
        if H_euler=={}:
            raise
    except:
        raise 

###################################################################
# GET THE EULER POLE VALUE
###################################################################
W=None

if args.plate is not None:
    if args.plate.strip() not in list( H_euler.keys() ):
        print("!!! ERROR: plate %s not found in euler file %s" % (args.plate,args.euler_file))
        sys.exit()
    else:
        W=value.reshape(3,1)

if args.value is not None:
    lpole = args.value.split('/')
    tmp_pole = lpole[-3:]
    W = np.array( list( map(float,tmp_pole) ) ).reshape(3,1)

print("-- Using Euler pole value: %8.4lf %8.4lf %8.4lf" % tuple(W.flatten().tolist()) )

###################################################################
# CASE IVEL
###################################################################

if args.ivel is not None:

    vel = VF.read(file_name=args.ivel)
    
    if verbose:
        print("-- Reading %s" % args.vel)
        print("-- Found %d in %s" % (vel.nsites(), args.vel))

    new_vel = vel.substract_pole(W, type_euler='euler')
    
    if args.ovel is not None:
        my_comment = ("Euler pole (%8.4lf %8.4lf %8.4lf) applied to %s" % tuple(W.flatten().tolist()+[args.ivel]) )
        import os
        if os.path.exists(args.ovel):
            print("!!! WARNING: %s exists. Results will be appended to it." % args.ovel)
        new_vel.write(args.ovel , comment= my_comment)
    else:
        for code in new_vel.lcode():
            new_vel.print_info_site(code, verbose=False)

###################################################################
# CASE PLL
###################################################################

if args.pll is not None:

    from pyacs.lib.gmtpoint import GMT_Point
    W[2,0] = -W[2,0] 
    vel = VF()

    # decipher pll option
    i=1
    for pll in args.pll:
        [lon,lat] = list( map(float, pll.split('/')[-2:]) )
        M=GMT_Point(code=("X%03d" % i), lon=lon,lat=lat,Ve=0.,Vn=0.,SVe=0.,SVn=0.,SVen=0.)
        vel.add_point(M)
        i = i + 1
    
    # prediction

    new_vel = vel.substract_pole(W, type_euler='euler')
    
    # print results
    if args.ovel is not None:
        my_comment = ("# Euler pole (%8.4lf %8.4lf %8.4lf) prediction" % tuple(W.flatten().tolist()) )
        new_vel.write(args.ovel , comment='# ')
    else:
        for code in new_vel.lcode():
            new_vel.print_info_site(code, verbose=False)
