#!/usr/bin/env python
'''
###################################################################
# SCRIPT    : pyacs_gvel_comb.py
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

# PYACS

import pyacs.lib.glinalg    
import pyacs.lib.coordinates
from pyacs.vel_field import Velocity_Field as VF

# SET STARTING TIME
t0 = time()

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="pyacs_gvel_comb.py\n"
prog_info+="--------------------------------------------------------------------------------------------------------------------------------------\n"
prog_info+="Combines GNSS velocity fields\n"
prog_info+="--------------------------------------------------------------------------------------------------------------------------------------\n"

prog_epilog="J.-M. Nocquet (Geoazur-IRD-CNRS-OCA-UCA & IPGP) & R. Tissandier (IPGP) - March 2019"


#----------------------------------------------------------
# REQUIRED INPUT FOR ALL INVERSION SCHEME
#----------------------------------------------------------

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,formatter_class=RawTextHelpFormatter)

parser.add_argument('-i', action='store', dest='input',type=str,required=True,help='text file including the list of GMT psvelo files to be combined. The first will be used to defined the reference. A reweighting factor can be provided.')
parser.add_argument('-e', action='store',dest='experiment',type=str,required=True,help='experiment name')
parser.add_argument('--verbose', '-v', action='count', default=0, help='verbose mode')
parser.add_argument('--debug', action='count', default=0, help='debug mode')


###################################################################
# ARGS PARSING
###################################################################


def main():
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()

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

    print("-- STEP 1: Reading input velocity files")

    H_gvel = {}
    H_gvel_name = {}
    H_gvel_sigma = {}

    n_gvel = 0
    fi = open(args.input)

    for line in fi:
        line = line.strip()
        if line[0] == '#':
            print("-- COMMENT: ", line)
            continue

        lline = line.split()

        if len(lline) > 1:
            H_gvel_sigma[n_gvel] = float(lline[1])

        print("-- reading ", lline[0])
        try:
            H_gvel_name[n_gvel] = lline[0]
            H_gvel[n_gvel] = VF.read(file_name=lline[0], verbose=False)
        except Exception:
            print("!!!ERROR. Could not read ", lline[0])
            sys.exit()

        print("-- found %s sites " % H_gvel[n_gvel].nsites())

        if H_gvel[n_gvel].nsites() != len(H_gvel[n_gvel].lcode()):
            print("!!!ERROR: There are duplicated sites in ", lline[0])
            sys.exit()

        n_gvel = n_gvel + 1

    print("-- Number of velocity fields to be combined: ", n_gvel)

    print("-- STEP 2: Getting the sites to be estimated")

    lsite = []
    for gvel in H_gvel.values():
        lsite = lsite + gvel.lcode()
    n_gps = len(lsite)

    print("-- Number of observation: ", n_gps * 2)

    lsite = list(set(sorted(lsite)))
    np_lsite = np.array(lsite)

    n_comb_gps = len(lsite)
    print("-- Number of sites in combined solution: ", n_comb_gps)

    print("-- STEP 3: Building the linear system")

    G = np.zeros((2 * n_gps + 3, 2 * n_comb_gps + 3 * n_gvel))
    d_obs = np.zeros(2 * n_gps + 3)
    VCV_d_obs = np.zeros((2 * n_gps + 3, 2 * n_gps + 3))

    coor = np.zeros((n_comb_gps, 2))

    idx_obs = 0
    for idx_gvel, gvel in H_gvel.items():
        for M in gvel.sites:
            if debug:
                print("  -- processing site %s from solution %d " % (M.code, idx_gvel))

            R = pyacs.lib.coordinates.mat_rot_general_to_local(np.radians(M.lon), np.radians(M.lat))
            (x, y, z) = pyacs.lib.coordinates.geo2xyz(np.radians(M.lon), np.radians(M.lat), M.he)
            Ai = np.zeros([3, 3], float)
            Ai[0, 1] = z
            Ai[0, 2] = -y
            Ai[1, 0] = -z
            Ai[1, 2] = x
            Ai[2, 0] = y
            Ai[2, 1] = -x
            Ai = Ai / 1.E6

            RAi = np.dot(R, Ai)[:2, :]

            idx_M_in_G = np.where(np_lsite == M.code)[0][0]
            G[idx_obs*2:idx_obs*2+2, idx_M_in_G*2:idx_M_in_G*2+2] = np.eye(2, 2)
            G[idx_obs*2:idx_obs*2+2, n_comb_gps*2+idx_gvel*3:n_comb_gps*2+idx_gvel*3+3] = RAi

            coor[idx_M_in_G, :] = [M.lon, M.lat]

            d_obs[idx_obs*2] = M.Ve
            d_obs[idx_obs*2+1] = M.Vn

            CORR = np.eye(2, 2)
            CORR[0, 1] = M.SVen
            CORR[1, 0] = M.SVen
            STD = np.array([M.SVe, M.SVn])
            VCV = pyacs.lib.glinalg.corr_to_cov(CORR, STD)

            if idx_gvel in H_gvel_sigma.keys():
                VCV = VCV * H_gvel_sigma[idx_gvel]**2

            if np.min(np.diag(VCV)) < 0:
                sys.exit()

            VCV_d_obs[idx_obs*2, idx_obs*2] = VCV[0, 0]
            VCV_d_obs[idx_obs*2+1, idx_obs*2+1] = VCV[1, 1]
            VCV_d_obs[idx_obs*2, idx_obs*2+1] = VCV[0, 1]
            VCV_d_obs[idx_obs*2+1, idx_obs*2] = VCV[1, 0]

            idx_obs = idx_obs + 1

    rf_constraints = 1.E-3

    G[idx_obs*2:idx_obs*2+3, n_comb_gps*2:n_comb_gps*2+3] = np.eye(3, 3)
    VCV_d_obs[idx_obs*2, idx_obs*2] = rf_constraints**2
    VCV_d_obs[idx_obs*2+1, idx_obs*2+1] = rf_constraints**2
    VCV_d_obs[idx_obs*2+2, idx_obs*2+2] = rf_constraints**2

    print("-- STEP 4: Solving the linear system")

    SOL, COV_SOL, RESIDUALS, chi2 = pyacs.lib.glinalg.lscov_full(G, d_obs, VCV_d_obs, verbose=verbose)

    print(" -- Global chi2 of combination: %.2lf " % chi2)

    vel_comb_ve_vn = SOL[:n_comb_gps*2].reshape(-1, 2)
    vcv_vel_comb_ve_vn = COV_SOL[:n_comb_gps*2, :n_comb_gps*2]

    corr, sigma_m = pyacs.lib.glinalg.cov_to_corr(vcv_vel_comb_ve_vn)
    my_corr = np.append(np.diagonal(corr, offset=1), 0).reshape(-1, 2)[:, 0]

    COMB = np.append(vel_comb_ve_vn, sigma_m.reshape(-1, 2), axis=1)
    COMB = np.append(COMB, my_corr.reshape(-1, 1), axis=1)

    COOR = np.array(["%10.5f" % w for w in coor.reshape(coor.size)]).reshape(-1, 2)
    TMP = np.array(["%10.3f" % w for w in COMB.reshape(COMB.size)]).reshape(-1, 5)
    TMP = np.append(COOR, TMP, axis=1)
    COMB = np.append(TMP, np_lsite.reshape(-1, 1), axis=1)

    comb_file = args.experiment + '_comb.gmt'
    np.savetxt(comb_file, COMB, fmt='%10s %10s %10s %10s %10s %10s %10s %10s')

    W = SOL[n_comb_gps*2:].reshape(-1, 3)
    rot_file = args.experiment + '_rot.dat'
    np.savetxt(rot_file, W, fmt='%10.4f %10.4f %10.4f', header="rotation rate Wx,Wy,Wz of individual solutions. 1.E-9 rad/yr")

    tol_res = 0.01
    idx_obs = 0
    stat_file = args.experiment + '_stat.gmt'

    H_gvel_chi2 = {}
    H_gvel_wrms = {}
    H_gvel_nsite = {}
    H_gvel_ncsite = {}
    H_gvel_apvf = {}
    H_gvel_rfn = {}

    for idx_gvel, gvel in H_gvel.items():
        print("-- Residuals for solution: %d %s" % (idx_gvel, H_gvel_name[idx_gvel]))

        H_gvel_nsite[idx_gvel] = 0
        H_gvel_ncsite[idx_gvel] = 0
        chi2_sol = 0.
        n_res_sol = 0
        denom = 0.

        res_file = args.experiment + '_' + H_gvel_name[idx_gvel].split('/')[-1].split('.')[0] + '_res.dat'
        print("-- saving residuals in %s" % res_file)
        H_gvel_rfn[idx_gvel] = res_file

        rf = open(res_file, 'w')

        for M in gvel.sites:
            H_gvel_nsite[idx_gvel] = H_gvel_nsite[idx_gvel] + 1

            if np.sqrt(RESIDUALS[2*idx_obs]**2 + RESIDUALS[2*idx_obs+1]**2) > tol_res:
                H_gvel_ncsite[idx_gvel] = H_gvel_ncsite[idx_gvel] + 1
                rf.write("%10.5lf %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %s\n" % (
                    M.lon, M.lat, RESIDUALS[2*idx_obs], RESIDUALS[2*idx_obs+1], M.SVe, M.SVn, M.SVen, M.code))
                chi2_sol = chi2_sol + (RESIDUALS[2*idx_obs]**2/M.SVe**2 + RESIDUALS[2*idx_obs+1]**2/M.SVn**2)
                n_res_sol = n_res_sol + 1
                denom = denom + 1./M.SVe**2 + 1./M.SVn**2

            idx_obs = idx_obs + 1

        H_gvel_chi2[idx_gvel] = chi2_sol
        if n_res_sol != 0:
            H_gvel_wrms[idx_gvel] = np.sqrt(chi2_sol/(2*n_res_sol))
        else:
            H_gvel_wrms[idx_gvel] = 0.0
        if denom != 0:
            H_gvel_apvf[idx_gvel] = np.sqrt(chi2_sol/denom)
        else:
            H_gvel_apvf[idx_gvel] = 0.0
        rf.close()

        print("-- individual chi2: %.2lf " % chi2_sol)
        print("-- individual post. var. factor: %.2lf" % H_gvel_apvf[idx_gvel])
        print("-- individual wrms: %.2lf mm/yr " % H_gvel_wrms[idx_gvel])

    sf = open(stat_file, 'w')
    my_line = "----------------------------------------------------------------------------------------------------------\n"
    sf.write(my_line)
    sf.write("n                sol_name                  nsite    ncsite      wrms      sigma      psigma   sigma*psigma\n")
    sf.write(my_line)

    for idx_gvel, gvel in H_gvel.items():
        rfn = H_gvel_rfn[idx_gvel]
        nsite = H_gvel_nsite[idx_gvel]
        ncsite = H_gvel_ncsite[idx_gvel]
        wrms = H_gvel_wrms[idx_gvel]
        apvf = H_gvel_apvf[idx_gvel]
        sigma = H_gvel_sigma[idx_gvel]
        sf.write("%02d %-40s %04d     %04d %10.2lf %10.2lf %10.2lf %10.2lf\n" % (
            idx_gvel, rfn, nsite, ncsite, wrms, sigma, apvf, apvf*H_gvel_sigma[idx_gvel]))
    sf.write(my_line)
    sf.close()


if __name__ == "__main__":
    main()