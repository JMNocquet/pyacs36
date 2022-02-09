#!/usr/bin/env python

###################################################################
# SCRIPT    : pyacs_gvel_estimate_pole.py
# AUTHOR    : JM NOCQUET
# DATE      : September 2011 / October 2017
# INPUT     : gmt psvelo file
#           : experiment name
#           | input plate file
# OUTPUT    : 
# NOTE      :
#         
###################################################################

# Implemented linear optimization methods implemented 
# WLS = Weighted Least Squares
# LS = Least Squares (no weight)
# L1 = L1 Norm optimization 

l_opt_method = ('WLS', 'LS', 'L1')

###################################################################
# MODULES IMPORT
###################################################################

import argparse

from pyacs.lib.vel_field import Velocity_Field as VF

###################################################################
# PARSE ARGUMENT LINE
###################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-plate', action='store', type=str, dest='plate', required=True, help='plate file with format: euro = SJDV POST WTZR JOZE')
parser.add_argument('-vel', action='store', type=str, dest='vel', required=True, help='horizontal velocity field, GMT psvelo format')
parser.add_argument('-experiment', action='store', type=str, dest='experiment', required=True, help='experiment name')
parser.add_argument('-method', action='store', type=str, dest='method', default='L2', help="'L1' or 'L2' norm. default 'L2'")
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
# INIT LOG FILE
###################################################################

log_file = args.experiment + '.log'
main_log = open(log_file, 'w')

###################################################################
# READ GMT FILE
###################################################################

vel = VF.read(file_name=args.vel)

main_log.write("-- Reading %s\n" % args.vel)
main_log.write("-- Found %d sites in %s\n" % (vel.nsites(), args.vel))

if verbose:
    print("-- Reading %s" % args.vel)
    print("-- Found %d in %s" % (vel.nsites(), args.vel))

###################################################################
# READ INPUT POLE FILE
###################################################################

main_log.write("-- Reading input file for pole calculation: %s\n" % args.plate)
if verbose:
    print("-- Reading input file for pole calculation: %s" % args.plate)

H_PLATE = {}
fs = open(args.plate, 'r')
for line in fs:
    if (len(line) < 2):continue
    if (line[0] == '#'):continue
    lline = line.split('=')
    if len(lline) < 2:continue
    lplate = lline[0].split()
    lsite = lline[1]
    H_PLATE[lplate[0]] = lsite
fs.close()

main_log.write("-- Found: %d plate(s) defined in %s\n" % (len(H_PLATE.keys()) , args.plate))
if verbose:
    print("-- Found: %d plate(s) defined in %s " % (len(H_PLATE.keys()) , args.plate))

###################################################################
# LOOP ON PLATES
###################################################################

for plate in H_PLATE.keys():
    
    main_log.write("####################################################\n")
    main_log.write("-- Processing plate %s\n" % plate)
    
    print("\n-- Processing plate %s" % plate)
    
    linclude = H_PLATE[plate].split()
    
    main_log.write("-- # requested sites : %d\n" % len(linclude))
    main_log.write("-- list of requested sites : \n")
    nb_site_per_line = 10
    i = 0
    for site in linclude:
        main_log.write("%s " % site)
        i = i + 1
        if i > nb_site_per_line:
            main_log.write("\n" % linclude)
            i = 0
    main_log.write("\n" % linclude)

    plate_log = log_file + '.' + plate
    velplate = vel.subset(linclude)
    
    print("-- plate ", plate)
    # CHECK THAT AT LEAST 3 SITES CAN BE USED FOR CALCULATION        
    if (velplate.nsites() < 2):
        main_log.write("!!! ERROR : less than 2 sites requested for plate :%s \n", plate)
        main_log.write("!!! ERROR : Euler pole calculation not possible for plate : %s \n", plate)
        
        print("!!! ERROR : less than 2 sites requested for plate :%s", plate)
        print("!!! ERROR : Euler pole calculation not possible for plate : %s", plate)
        
        continue

    # CHECK THAT ALL ASKED SITES WERE FOUND        
    if (velplate.nsites() < len(linclude)):
        main_log.write("! WARNING : not all the requested sites were found for plate %s found in %s\n" % (plate, args.vel))
        print("! WARNING : not all the requested sites were found for plate %s found in %s" % (plate, args.vel))
        
        i = 0
        lmissing = [item for item in linclude if not item in velplate.lcode()]
        for site in lmissing:
            main_log.write("WARNING: missing site: %s \n" % site)
            i = i + 1
            if i > nb_site_per_line:
                main_log.write("\n" % linclude)
                i = 0
        main_log.write("\n" % linclude)
    
    # BUILD LINEAR SYSTEM
    main_log.write("-- Now doing Euler pole inversion\n")
    
    pole,vcv = velplate.pole(exp=plate, log=plate_log)
    
    # WRITING OUTPUTS
    
    # Euler pole residuals

    name_out = args.experiment + '_res_pole_' + plate
    main_log.write("-- writing results to: %s\n" % name_out)
    velplate.substract_pole(W=pole, type_euler='rot').write(name_out, verbose=verbose, comment='')
    
    # New velocity field wrt the new pole
    name_out = args.experiment + '_res_' + plate
    main_log.write("-- writing results to: %s\n" % name_out)
    vel.substract_pole(W=pole,type_euler='rot').write(name_out, verbose=verbose, comment='')
main_log.close()

