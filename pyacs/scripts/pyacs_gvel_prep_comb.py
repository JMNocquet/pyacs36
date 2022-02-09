#!/usr/bin/env python
'''
###################################################################
# SCRIPT    : pyacs_prep_gvel_comb.py
# AUTHOR    : 
# DATE      : March 2019
# INPUT     :  
# OUTPUT    : 
# NOTE      :
###################################################################
'''


###################################################################
# MODULES IMPORT
###################################################################

# GENERAL

import sys, os
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
from time import time
import linecache 

# PYACS

import pyacs.lib.glinalg    
import pyacs.lib.coordinates
from pyacs.lib.vel_field import Velocity_Field as VF

# SET STARTING TIME
t0 = time()

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="pyacs_prep_gvel_comb.py\n"
prog_info+="--------------------------------------------------------------------------------------------------------------------------------------\n"
prog_info+="Prepares the input files for the combination of GNSS velocity fields\n"
prog_info+="--------------------------------------------------------------------------------------------------------------------------------------\n"

prog_epilog="J.-M. Nocquet (Geoazur-IRD-CNRS-OCA-UCA & IPGP) & R. Tissandier (IPGP) - March 2019"


#----------------------------------------------------------
# REQUIRED INPUT FOR THE PREPARATION
#----------------------------------------------------------

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,formatter_class=RawTextHelpFormatter)

parser.add_argument('-i', action='store', dest='input',type=str,required=True,help='text file including the list of GMT psvelo files to be combined. The first one has to be the reference.')
parser.add_argument('-d', action='store',dest='critic_distance',type=float,required=True,help='critical distance in km')
parser.add_argument('--verbose', '-v', action='count', default=0, help='verbose mode')
parser.add_argument('--debug', action='count', default=0, help='debug mode')


if (len(sys.argv)<2):parser.print_help();sys.exit()

###################################################################
# ARGS PARSING
###################################################################

args = parser.parse_args()

# verbose
if args.verbose>0:
    print("-- Verbose mode")
    verbose=True
else:
    verbose=False
# debug
if args.debug>0:
    print("-- Debug mode")
    verbose=True
    debug=True
else:
    debug=False


###################################################################
# READ INPUT VELOCITY FIELDS
###################################################################

print("-- STEP 1: Reading the first input velocity file")

H_gvel_ref = {}
H_gvel = {} # dictionary of numpy 2D arrays including the velocity field (GMT psvelo format) without site name
H_gvel_site = {} # dictionary of numpy 1D arrays including the site name
H_gvel_file ={}
n_gvel = 0 # number of gvel to compare
rename_idx = 0 # index for renaming will be X000, X001 etc

fi = open(args.input) 

for line in fi:
    # strip
    line = line.rstrip() ##### strip=returns a list of all the words in the string
    # if commented then forget it
    if line[0] == '#':
        print("-- COMMENT: " , line)
        next()

    # reads line
    vel_name = line #### ligne du list txt
    ##### split=returns a list of all the words in the string

    # tries to read the velocity field
    print("-- reading: %s" % vel_name)
    try:
        H_gvel_file[ n_gvel ] = vel_name
        H_gvel[ n_gvel ] = np.genfromtxt(vel_name,comments='#',usecols=(0,1,2,3,4,5,6))
        H_gvel_site[ n_gvel ] = np.genfromtxt(vel_name,comments='#',usecols=(7) , dtype=str)
    except:
        print("!!!ERROR. Could not read " , vel_name)
        sys.exit()

    print("-- found %s sites " % len(H_gvel_site[ n_gvel ]) )

    # check that all site name are uniq
    while np.unique( H_gvel_site[ n_gvel ] ).size != H_gvel_site[ n_gvel ].size:
        print("!!! There are duplicated codes in %s" % H_gvel_file[ n_gvel ])
        lduplicated=list( set([x for x in H_gvel_site[n_gvel].tolist() if H_gvel_site[n_gvel].tolist().count(x) > 1]) )
        for site in sorted( lduplicated ):
            idx = np.argwhere(H_gvel_site[ n_gvel ]==site)[-1,0]
            new_site = ("X%03d" % rename_idx)
            print("!!! %s: renaming duplicated site %s as %s" % (H_gvel_file[n_gvel],site,new_site) )
            H_gvel_site[n_gvel][idx] = new_site
            rename_idx = rename_idx + 1

    # dealing with close sites within the same solution
    # there is no obvious solution.
    # mean solution

    lsite_remove=[]
    
    for si in np.arange( H_gvel_site[ n_gvel ].size ):
        for sj in np.arange( si+1, H_gvel_site[ n_gvel ].size ):
            coo_i = H_gvel[n_gvel][si,:]
            coo_j = H_gvel[n_gvel][sj,:]

            diff_coo = coo_j[:2] - coo_i[:2] 
            diff_km = np.sqrt(np.sum(diff_coo**2))*111.
    
            if diff_km <= args.critic_distance and H_gvel_site[n_gvel][si] != H_gvel_site[n_gvel][sj]:
                print("!!! site %s in solution %s and site %s in %s are %.2lf km apart (<%.2lf):" % \
                      (H_gvel_site[n_gvel][si],H_gvel_file[n_gvel],H_gvel_site[n_gvel][sj],H_gvel_file[n_gvel],diff_km,args.critic_distance))
                print("!!! %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf" % tuple(coo_i))
                print("!!! %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf" % tuple(coo_j))
                print("!!! Taking the mean position and value")
                ve = ( coo_i[2]/coo_i[4]**2 + coo_j[2]/coo_j[4]**2) / ( 1./coo_i[4]**2 + 1./coo_j[4]**2) 
                vn = ( coo_i[3]/coo_i[5]**2 + coo_j[3]/coo_j[5]**2) / ( 1./coo_i[5]**2 + 1./coo_j[5]**2)
                coo_i[2] = ve
                coo_i[3] = vn
                H_gvel[n_gvel][si,:] = coo_i
                print("!!! %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf" % tuple(coo_i))
                print("!!! %s will be removed from solution: %s" % (H_gvel_site[n_gvel][sj],H_gvel_file[n_gvel]))
                lsite_remove.append( H_gvel_site[n_gvel][sj] )
    
    # remove sites
    if lsite_remove != []:
        print("!!! %d sites will be removed from solution: %s" % (len(lsite_remove), H_gvel_file[n_gvel] ))
        print("!!! " , lsite_remove)
        lidx_remove=[]
        for rsite in lsite_remove:
            lidx_remove.append( np.argwhere( H_gvel_site[ n_gvel ]== rsite )[-1,0] )
    
        H_gvel_site[ n_gvel ] = np.delete( H_gvel_site[ n_gvel ] , lidx_remove )
        H_gvel[ n_gvel ]      = np.delete( H_gvel [ n_gvel ] , lidx_remove , axis=0 )
    
        print("-- keeping %s sites " % len(H_gvel_site[ n_gvel ]) )
    
    n_gvel = n_gvel + 1


print("-- Number of velocity fields to be combined: " , n_gvel)



###################################################################
# TEST ON SITE NAME
###################################################################   


print("-- STEP 2: Checking site name and distances in solutions")

for i in np.arange( n_gvel ):
            
    for j in np.arange( i+1 , n_gvel ):
        print("-- checking files: %s vs %s" % ( H_gvel_file[i] , H_gvel_file[j]))

        # check that common name are actually the same sites
        
        np_common = np.intersect1d(H_gvel_site[i], H_gvel_site[j] )
        print("-- Found %d sites with the same name " % np_common.size)
        
        for s in np.arange( np_common.size ):
            site = np_common[s]
            idxi = np.argwhere( H_gvel_site[i]==site)[0,0]
            coo_i = H_gvel[i][idxi,:]

            idxj = np.argwhere( H_gvel_site[j]==site)[0,0]
            coo_j = H_gvel[j][idxj,:]

            #print(site , coo_j-coo_i)
            
            diff_coo = coo_j[:2] - coo_i[:2] 
            diff_km = (np.sqrt(np.sum(diff_coo**2))*111.)
            
            if diff_km > args.critic_distance:
                print("!!! Pb. site %s has the same name in solutions %s and %s but its coordinates are:" % (site,H_gvel_file[i],H_gvel_file[j]))
                print("!!! %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf" % tuple(coo_i))
                print("!!! %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf" % tuple(coo_j))

                new_site = ("X%03d" % rename_idx)
                print("!!! %s: not a common site renaming site %s as %s" % (H_gvel_file[j],site,new_site) )
                H_gvel_site[j][idxj] = new_site
                rename_idx = rename_idx + 1
                
        # results

        np_common = np.intersect1d(H_gvel_site[i], H_gvel_site[j] )
        print("-- After correction, found %d sites with the same name " % np_common.size)

        # check all distances

        for si in np.arange( H_gvel_site[i].size ):
            for sj in np.arange( H_gvel_site[j].size ):
                coo_i = H_gvel[i][si,:]
                coo_j = H_gvel[j][sj,:]

                diff_coo = coo_j[:2] - coo_i[:2] 
                diff_km = np.sqrt(np.sum(diff_coo**2))*111.
        
                if diff_km <= args.critic_distance and H_gvel_site[i][si] != H_gvel_site[j][sj]:
                    print("!!! site %s in solution %s and site %s in %s are %.2lf km apart (<%.2lf):" % \
                          (H_gvel_site[i][si],H_gvel_file[i],H_gvel_site[j][sj],H_gvel_file[j],diff_km,args.critic_distance))
                    print("!!! %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf" % tuple(coo_i))
                    print("!!! %10.5lf %10.5lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf" % tuple(coo_j))
    
                    print("!!! %s: close site renaming site %s as %s" % (H_gvel_file[j],H_gvel_site[j][sj],H_gvel_site[i][si]) )
                    H_gvel_site[j][sj] = H_gvel_site[i][si]
        
# save the results

for i in np.arange( n_gvel ):
    
    cln_file_name = H_gvel_file[i]+'.cln'
    
    print("-- Saving clean file: %s" % cln_file_name )
    
    f = open( cln_file_name , 'w')
    
    for si in np.arange( H_gvel_site[i].size ):
        f.write("%10.5lf %10.5lf  %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % tuple( H_gvel[i][si,:].tolist()+[H_gvel_site[i][si]]) )
        
    f.close()
        
        

# # crée un tableau de nombres pour les noms de stations à modifier
# 
# number = [a for a in np.arange(1000,1500,1)]
# n = 0
# 
# # trie le fichier de la plus petite à la plus grande longitude
# 
# sorted(H_gvel_ref, key = lambda colonnes: colonnes[0]) 
# sorted(H_gvel, key = lambda colonnes: colonnes[0])
# 
# print("Looking for common sites between solutions ")
# 
# # check for common names between the different files
# 
# for i in H_gvel_ref[ 7 ].names_ref():   
#     for j in H_gvel[ 7 ].names:   
#         if H_gvel_ref[ i ].name_ref() == H_gvel[ j ].name():
# 
#             # spherical distance between the common sites
# 
#             spherical_distance = pyacs.lib.coordinates.geo_spherical_distance(H_gvel_ref[ i ].lon_ref(),H_gvel_ref[ i ].lat_ref(),6.371E6,H_gvel[ j ].lon(),H_gvel[ j ].lon(),6.371E6,unit='dec_deg',Rt=6.371E6)
#            
#             # if sites separated by more than 1 km, it must be different sites
#             if spherical_distance >= critic_distance:
#                 H_gvel[ j ] = number[n] 
#                 n += 1
#             else:
#                 sys.exit()
#         
#         else:
#             #if le point se rapproche til d'un autre point existant ?
#                 #on lui met le nom de ce point existant H_gvel[ j ].names() = ??
#             else 
#             sys.exit()
#            
# 
# ####### MEMO ######
# # OBJETS = H_GVEL ET H_GVEL_REF contenant les VF soit 7 colonnes différentes ?
# # ATTRIBUTS =   lsite ?
# #               nsite ? 
# 
# 
# ###################################################################
# # PRINT RESULTS
# ###################################################################
# 
# #comment incrémenter tout ça dans le output ?
# 
# cleaned_file = args.output+'_cleaned.gmt'
# np.savetxt(cleaned_file, bla, fmt='%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f', header="cleaned file lon lat vE vN sE sN CorrelationEN site")
# 
# 

