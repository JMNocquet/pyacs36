#!/usr/bin/env python

###################################################################
# SCRIPT    : pyacs_gvel_strain.py
# AUTHOR    : JM NOCQUET
# DATE      : September 2011 / January 2018 / Python3 May 2019
# INPUT     : gmt psvelo file
#           : experiment name
#           : input strain file
# OUTPUT    : 
# NOTE      :
#         
###################################################################

# Implemented linear optimization methods implemented 
# WLS = Weighted Least Squares
# LS = Least Squares (no weight)
# L1 = L1 Norm optimization 

l_opt_method=('WLS')


###################################################################
# MODULES IMPORT
###################################################################

import argparse, sys, os
from os import path,mkdir

from pyacs.vel_field import Velocity_Field as VF

###################################################################
# PARSE ARGUMENT LINE
###################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--strain', action='store', type=str, dest='strain',default=None,help='strain file with format SJDV POST WTZR JOZE (1 list per line)')
parser.add_argument('--site', action='append', type=str, default=[],dest='lsite',help='sites to be used for strain rate calculation. Must be repeated at least 3 times')
parser.add_argument('-vel', action='store', type=str, dest='vel',required=True,help='horizontal velocity field, GMT psvelo format')
parser.add_argument('-experiment', action='store', type=str, dest='experiment',required=True,help='experiment name')
parser.add_argument('-method', action='store', type=str, dest='method',default='L2',help="'L1' or 'L2' norm. default 'L2'" )
parser.add_argument('--verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('--debug', action='count',default=0,help='debug mode')


def main():
    args = parser.parse_args()

    if args.verbose > 0:
        print("-- Verbose mode")
        verbose = True
    else:
        verbose = False

    if args.debug > 0:
        print("-- Debug mode")
        verbose = True
        debug = True
    else:
        debug = False

    if (args.strain) is None and (args.lsite == []):
        print('!!! ERROR: either --site or --strain option are required.')
        sys.exit()

    log_file = args.experiment + '.log'
    main_log = open(log_file, 'w')

    vel = VF.read(file_name=args.vel)

    main_log.write("-- Reading %s\n" % args.vel)
    main_log.write("-- Found %d sites in %s\n" % (vel.nsites(), args.vel))

    if verbose:
        print("-- Reading %s" % args.vel)
        print("-- Found %d in %s" % (vel.nsites(), args.vel))

    i = 0
    H_STRAIN = {}

    if (args.strain) is not None:
        main_log.write("-- Reading input file for strain calculation: %s\n" % args.strain)
        if verbose:
            print("-- Reading input file for pole calculation: %s" % args.strain)

        fs = open(args.strain, 'r')
        for line in fs:
            if (len(line) < 2):
                continue
            if (line[0] == '#'):
                continue
            lsite = line.split()
            if len(lsite) < 3:
                print('! WARNING: a minimum of three sites is required: ', line)
                main_log.write('! WARNING: a minimum of three sites is required: %s\n' % line)
            else:
                i = i + 1
                H_STRAIN[i] = lsite
        fs.close()

        main_log.write("-- Found: %d strain rate to be estimated in %s\n" % (i, args.strain))
        if verbose:
            print("-- Found: %d strain rate defined in %s " % (i, args.strain))

    if args.lsite != []:
        if len(args.lsite) < 3:
            print('! WARNING: a minimum of three sites is required: ', ' '.join(args.lsite))
            main_log.write('! WARNING: a minimum of three sites is required: %s\n' % (' '.join(args.lsite)))
            print('! WARNING: skipping this subset for strain rate calculation')
            main_log.write('! WARNING: skipping this subset for strain rate calculation\n')
        else:
            i = i + 1
            H_STRAIN[i] = args.lsite

    for index in list(H_STRAIN.keys()):
        main_log.write("####################################################\n")
        main_log.write("-- Processing strain %d / %d\n" % (index, len(list(H_STRAIN.keys()))))

        print("-- Processing strain %d" % index)

        linclude = H_STRAIN[index]

        main_log.write("-- # requested sites : %d\n" % len(linclude))
        main_log.write("-- list of requested sites : \n")
        nb_site_per_line = 10
        i = 0

        for site in linclude:
            main_log.write(" %s " % site)
            i = i + 1
            if i > nb_site_per_line:
                main_log.write("\n")
                i = 0
        main_log.write("\n")

        strain_log = ('strain_%04d.dat' % index)

        velstrain = vel.subset(linclude)

        if (velstrain.nsites() < len(linclude)):
            main_log.write("! WARNING : not all the requested sites were found for strain %i found in %s\n" % (index, args.vel))
            print("! WARNING : not all the requested sites were found for strain %i found in %s" % (index, args.vel))
            i = 0
            lmissing = [item for item in linclude if item not in velstrain.lcode()]
            for site in lmissing:
                main_log.write("WARNING: missing site: %s \n" % site)
                i = i + 1
                if i > nb_site_per_line:
                    main_log.write("\n")
                    i = 0
            main_log.write("\n")

        main_log.write("-- Now doing strain rate inversion\n")

        velstrain.strain(velstrain.lcode(), save=strain_log, method='WLS', verbose=verbose)

        main_log.write("-- writing results to: %s\n" % strain_log)
        print("-- writing results to: %s" % strain_log)

    print('-- log written in ', log_file)
    main_log.close()


if __name__ == "__main__":
    main()
