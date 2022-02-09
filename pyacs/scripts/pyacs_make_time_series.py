#!/usr/bin/env python

###################################################################
# SCRIPT    : pyacs_time_series.py
# AUTHOR    : JM NOCQUET
# DATE      : February 2011
# INPUT     : list of sinex files
#           : reference sinex
#           : experiment name
#           | options file
# OUTPUT    : 
# NOTE      :
#         
###################################################################
# TODO
# file Helmert_outliers and Helmert_summary
# ref_bias_pos & ref_bias_vel

###################################################################
# MODULES IMPORT
###################################################################


import time,os,sys,datetime
import numpy as np
import argparse
import shutil
import glob

from colors import red

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


import pyacs.lib
import pyacs.sol.log
import pyacs.sol.sinex as SSinex
import pyacs.sol.helmert as Helmert
import pyacs.sol.discontinuity as Discontinuities
import pyacs.sol.read_conf as Read_Conf_Pyacs

from pyacs.sinex.sinex import sinex
from pyacs.sinex import snxutils

from pyacs.gts.Gts import Gts
from pyacs.gts.Sgts import Sgts

###################################################################
# routine H_common from pck
###################################################################
def get_H_common( ref_ts, H_Gpoint_free, prefit=3.5):
    # import
    from pyacs.sol.gpoint import Gpoint
    import pyacs.lib.astrotime as at
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # initialize
    H_common_Gpoint = {}
    Hdistance = {}
    strict = True

    # get free sinex epoch in datetime using the first point
    [code,soln] = list(H_Gpoint_free.keys())[0]

    datetime_freesinex = at.decyear2datetime( H_Gpoint_free[code,soln].epoch )
    str_datetime_freesinex = datetime_freesinex.isoformat()[:10]
    DEBUG( datetime_freesinex )
    # find common sites for the given date
    lcode_ref = sorted(ref_ts.lcode() )
    DEBUG(lcode_ref)
    # loop
    for code, soln in list(H_Gpoint_free.keys()):
        if code in lcode_ref:
            DEBUG("%s %s " % (code,soln))
            M = H_Gpoint_free[code, soln]
            DEBUG("free M.X: %.10lf " % M.X )
            # convert ts to pandas time series, this can be done earlier once for performance improvement
            pd_ts = ref_ts.__dict__[code]\
            #.to_pandas_df( data_xyz=True , uncertainty=True )

            # get XYZ position from ref time series
            try:
                [X_ref, Y_ref, Z_ref, SX_ref, SY_ref, SZ_ref] = pd_ts.loc[str_datetime_freesinex].to_numpy().flatten()[:6]
            except:
            #if pd_ts.loc[ str_datetime_freesinex ].empty:
                WARNING("No data for date %s for site %s " %(str_datetime_freesinex,code))
                #DEBUG(pd_ts.loc[ str_datetime_freesinex ])
                #DEBUG(pd_ts.loc[ str_datetime_freesinex ].to_numpy())
                #DEBUG("pck M.X: %15.8lf " % X_ref )
                continue

            N = Gpoint(X=X_ref, Y=Y_ref, Z=Z_ref, \
                   SX=SX_ref, SY=SY_ref, SZ=SZ_ref, \
                   VX=0., VY=0., VZ=0., \
                   SVX=0., SVY=0., SVZ=0., \
                   epoch=M.epoch, code=code, pt='A', soln=soln)

            distance = M.xyz_distance(N)
            if (distance > prefit):
                WARNING("Bad prefit for %s %10.1lf meters (threshold %8.1lf) " % (M.code, M.xyz_distance(N), prefit ))

            else:
                Hdistance[code, soln] = distance
                H_common_Gpoint[code, soln] = N

    if strict:
        median_distance = np.median(list(Hdistance.values()))

        VERBOSE(("Median prefit: %10.2lf mm" % (median_distance * 1.E3)))

        for code, soln in list(Hdistance.keys()):
            if Hdistance[code, soln] > 50. * median_distance:
                VERBOSE("Deleting %s %s bad prefit=%.3lf m" %( code, soln, Hdistance[code, soln] * 1.E3))
                del H_common_Gpoint[code, soln]

    return (H_common_Gpoint)


###################################################################
# START TIME
###################################################################
start = time.time()   
str_date = datetime.datetime.now().strftime("%Y_%m_%d")

###################################################################
# PARSE ARGUMENT LINE
###################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-lsinex', action='store', type=str, dest='lsinex_name',required=True,help='file including the list of sinex/glx files to be processed')
parser.add_argument('-experiment', action='store', type=str, dest='expt',required=True,help='experiment name')
parser.add_argument('-ref_sinex', action='store', type=str, dest='ref_sinex_name',default=None,help='reference sinex used to define the reference frame')
parser.add_argument('--ref_apr', action='store', type=str, dest='ref_apr_name',default=None,help='reference apr file used to define the reference frame')
parser.add_argument('--ref_pck', action='store', type=str, dest='ref_pck',default=None,help='reference pck time series used to define the reference frame')
parser.add_argument('--eq_rename', action='store', type=str, dest='eq_rename',default=None,help='globk eq_rename file. Will be used with glred if glx are provided')
parser.add_argument('--discontinuity', action='store', type=str, dest='discontinuity',default=None,help='discontinuity file (IGS format ftp://igs-rf.ign.fr/pub/discontinuities/soln.snx)')
parser.add_argument('--codomes', action='store', type=str, dest='codomes',default=None,help='codomes file (IGS format ftp://igs-rf.ign.fr/pub/DOMES/codomes_gps_coord.snx)')
parser.add_argument('--psd', action='store', type=str, dest='psd',default=None,help='post-seismic deformation file (IGS format ftp://igs-rf.ign.fr/pub/psd/psd_IGS.snx)')
parser.add_argument('--conf_file', action='store', type=str, dest='options_file', default=None, help='configuration file for pyacs')
parser.add_argument('--dikin_outlier', action='store', type=float, dest='dikin_outlier', default=3, help='threshold value for outliers rejection after L1 estimation (1: very strict, 3: standard, 10: loose)')
parser.add_argument('--rep', action='store', type=float, dest='min_repeatability', default=None, help='relative weight given to the minimum repeatability condition; default 0.')
parser.add_argument('--method', action='store', type=str, dest='method', default='Dikin_LS', help='Estimation method Dikin_LS or LS')
parser.add_argument('--verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('--debug', action='count',default=0,help='debug mode')
parser.add_argument('--replace', action='count',default=0,help='replace mode')
parser.add_argument('--uncertainty', action='count',default=0,help='will calculate daily coordinates uncertainties. Requires either sinex or glx files are provided in lsinex.')
parser.add_argument('--glred', action='count',default=0,help='will run glred together with pyacs to get the time series. Requires glx files are provided in lsinex.')
parser.add_argument('--no_clean', action='count',default=0,help='do not remove created sinex files.')
parser.add_argument('--pck_only', action='count',default=0,help='only writes pck file.')


args = parser.parse_args()

###################################################################
# OPTIONAL ARGUMENTS
###################################################################

# verbose
if args.verbose>0:
    MESSAGE("Verbose mode")
    verbose=True
else:
    verbose=False

# debug
if args.debug>0:
    MESSAGE("Debug mode")
    verbose=True
    debug=True
else:
    debug=False

# replace
if args.replace>0:
    MESSAGE("Replace mode")
    replace=True
else:
    replace=False

# uncertainty

if args.uncertainty>0:
    MESSAGE("Uncertainty mode")
    uncertainty=True
else:
    uncertainty=False

# remove sinex files
if args.no_clean>0:
    MESSAGE("Will keep sinex files")
    clean = False
else:
    clean = True
    
# glred option
if args.glred:
    MESSAGE("Will run glred")
    glred_opt = True
else:
    glred_opt = False

# pck only
if args.pck_only>0:
    MESSAGE("pck only option. Will only writes a pck file.")
    pck_only=True
else:
    pck_only = False

# save command line
cmd_line = " ".join(sys.argv)

###################################################################
# SET VERBOSE MODE THROUGH LOGGING MODULE
###################################################################

logging.getLogger("my_logger").setLevel(logging.WARNING)

if not(verbose):
    # WARNING LEVEL: ONLY MAJOR STEPS AND WARNINGS WILL BE PRINTED
    logging.getLogger("my_logger").setLevel(logging.WARNING)
else:
    if verbose:
        # VERBOSE MODE
        logging.getLogger("my_logger").setLevel(logging.INFO)
    else:
        # WARNING MODE
        logging.getLogger("my_logger").setLevel(logging.WARNING)

if debug:
    # DEBUG MODE
    logging.getLogger("my_logger").setLevel(logging.DEBUG)

if logging.getLogger("my_logger").getEffectiveLevel() == logging.WARNING:
    MESSAGE("verbose level: WARNING ONLY")

if logging.getLogger("my_logger").getEffectiveLevel() == logging.INFO:
    MESSAGE("verbose level: VERBOSE")

if logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG:
    MESSAGE("verbose level: DEBUG")


###################################################################
# READS OPTIONS FILE
###################################################################
if args.options_file:
    conf=Read_Conf_Pyacs.Conf(args.options_file,verbose=verbose)
    DEBUG("conf.rename: %s " % conf.rename)
else:
    conf=Read_Conf_Pyacs.Conf()

if conf.min_repeat != []:
    previous_sinex=SSinex.SSinex()

###################################################################
# READS FILE INCLUDING THE SINEX FILES NAME TO BE PROCESSED
###################################################################

MESSAGE('Reading list of files to be processed: %s' % args.lsinex_name)

# check lsinex has some files to process
try:
    lsinex=np.genfromtxt(args.lsinex_name,comments='#',dtype=str).tolist()
    if lsinex == []:
        from pyacs.sol.errors import PyacsSol_BadLsnx
        raise PyacsSol_BadLsnx( args.lsinex_name )
except PyacsSol_BadLsnx as error:
    # print PYACS ERROR
    ERROR( error, exit=True )

# ensure each file line is uniq
if np.size(lsinex)>1:
    lsinex=sorted(list(set(lsinex)))
else:
    lsinex=[lsinex]


# type_sol glx, snx or ssc

type_sol = 'unknown'

ll=lsinex[0].split('/')
sinex_basename=ll[-1]

suffix=sinex_basename.split('.')[-1]
if suffix in ['glx','GLX']:
    type_sol = 'glx'

if suffix in ['snx','SNX']:
    type_sol = 'snx'

if suffix in ['ss','SS','ssc','SSC']:
    type_sol = 'ssc'

###################################################################
# CHECK OPTIONS CONSISTENCY
###################################################################

if args.uncertainty and ( type_sol == 'ssc' ):
    ERROR("ERROR in arguments: uncertainty option requires either glx or sinex solutions as input",exit=True)
    sys.exit()

if glred_opt and ( type_sol != 'glx' ):
    ERROR("ERROR in arguments: glred option requires glx solutions as input")
    sys.exit()

MESSAGE("Number of solutions to be processed: %d" % (len(lsinex) ))
MESSAGE("Solution type: %s" % (type_sol ))

###################################################################
# OUTPUT DIRECTORY
###################################################################

wdir = args.expt

if not replace:

    ###################################################################
    # CREATE OUTPUT DIRECTORY
    ###################################################################
    
    VERBOSE("Creating working directory: %s " % args.expt)

    if os.path.isdir(wdir):
        ERROR(("%s directory already exists" % wdir),exit=True )
        sys.exit()
    
    elif os.path.isfile(wdir):
        ERROR(("%s file already exists" % wdir),exit=True )
        sys.exit()


    else:
        (path_stat,path_res_pos,path_maps,path_pos, path_pck, path_txyz,path_tsg,path_prt,path_lgred_pos,path_log,path_sinex,path_ssc) = \
        pyacs.sol.log.output_directories(wdir)
    

else:

    VERBOSE("Reading working directory: %s " % args.expt)
    
    (path_stat,path_res_pos,path_maps,path_pos, path_pck, path_txyz,path_tsg,path_prt,path_lgred_pos,path_log,path_sinex,path_ssc) = \
        pyacs.sol.log.output_directories(wdir)


# Change JMN 03/09/2019. There is no need to delete them
# if verbose: print("-- Cleaning pos and glred_pos directory of ", args.expt)
#    try:
#        shutil.rmtree(path_pos)
#        shutil.rmtree(path_pck)
#        shutil.rmtree(path_lgred_pos)
#    except:
#        pass
#    
#    os.mkdir(path_pos)
#    os.mkdir(path_lgred_pos)


###################################################################
# TSG OUTPUT FILES
###################################################################

if replace:
    sum_helmert_file=("%s/residuals_tsg_%s.dat" % (path_stat,str_date))
    tsg_file=("%s/tsg_%s.dat" % (path_tsg,str_date))
else:
    sum_helmert_file=path_stat+'/residuals_tsg.dat'
    tsg_file = path_tsg+'tsg.dat'

###################################################################
# VARIABLES INITIALIZATION
###################################################################

# TSG list as 2D numpy array
TSG=np.zeros((len(lsinex),8))

# list of sites used to define the frame
l_ref_code=[]

# dictionary of Helmert estimation statistics
H_tsg_stat={}
        

###################################################################
# READS REFERENCE SINEX FILE
###################################################################

if args.ref_sinex_name is not None:

    # USING PYACS ROUTINE
    MESSAGE("Reading reference SINEX file: %s " % args.ref_sinex_name)
    
    ref_sinex=SSinex.SSinex(name=args.ref_sinex_name)
    ref_sinex.read_section_estimate(lexclude=conf.ref_exclude_all , lonly=conf.ref_only_all)
    
    MESSAGE("Found %d code/soln in %s " % (len(list(ref_sinex.estimates.keys())), ref_sinex.name))
    
    # USING P. REBISCHUNG SINEX LIBRARY
#    print '-- Reading reference SINEX file using P. Rebischung sinex library: ',args.ref_sinex_name
#    ref_sinex = sinex.read(args.ref_sinex_name,dont_read =['matrices','comments','apriori','metadata'])

#    ref_sinex=Sinex.Sinex(args.ref_sinex_name)
#    ref_sinex.read_section_estimate(lexclude=conf.ref_exclude_all)
    
if args.ref_apr_name:

    MESSAGE("Reading reference apr file: %s"% args.ref_apr_name)
    ref_sinex=SSinex.SSinex(args.ref_apr_name)
    ref_sinex.read_apr(lexclude=conf.ref_exclude_all)

if args.ref_pck:
    MESSAGE("Reading reference time series pck file %s" % args.ref_pck )
    import pickle
    try:
        with open( args.ref_pck, "rb") as f:
            ref_pck = pickle.load(f)
        f.close()
    except:
        ERROR(("reading pck file: %s " % args.ref_pck),exit=True )
    MESSAGE("Converting time series from pyacs format to pandas time series")
    for code in ref_pck.lcode():
        ref_pck.__dict__[code] = ref_pck.__dict__[code].to_pandas_df(data_xyz=True, uncertainty=True)

    MESSAGE("Number of GNSS time series used as reference: %d " % ref_pck.n())

##################################################################################
# READS IGS DISCONTINUITY FILE - SINEX & DISCONTINUITY MUST BE MUTUALLY CONSISTENT
##################################################################################

discontinuities=Discontinuities.Discontinuities()
if args.discontinuity is not None:
    # USING PYACS ROUTINE
    MESSAGE("Reading IGS discontinuity file: %s " % args.discontinuity)
    discontinuities.read_igs_discontinuity(args.discontinuity)

    # USING P. REBISCHUNG SINEX LIBRARY
    #print '-- Reading IGS discontinuity file using P. Rebischung sinex library : ',args.discontinuity
    #solns = snxutils.read_solns('soln.snx')

##################################################################################
# CODOMES
##################################################################################

if args.codomes is not None:

    # USING P. REBISCHUNG SINEX LIBRARY
    MESSAGE("Reading IGS codomes using P. Rebischung sinex library %s "  % args.codomes)
    domes = snxutils.read_domes(args.codomes)

##################################################################################
# POST-SEISMIC DEFORMATION FILE IGS FORMAT PSD_IGS.snx
##################################################################################

if args.psd is not None:

    # USING P. REBISCHUNG SINEX LIBRARY
    MESSAGE("Reading IGS psd file using P. Rebischung sinex library: %s"  % args.psd)
    psd = sinex.read(args.psd)
else:
    psd = None


###################################################################
# STARTS LOOP FOR EACH SINEX FILES
###################################################################

# H_res_helmert
H_res_helmert = {}

# Sgts instance for reference sites residual time series
l_res_ts = Sgts(read=False)

# initialize sts
sts = Sgts(read=False)

##########################################################################
MESSAGE("START OF PROCESSING - LOOP ON SINEX FILES" , level=2)
##########################################################################

index_sinex = 0

for sinex_name in lsinex:
    
    # Helmert and sites covariance initialization
    T=None
    H_COV = None 

    MESSAGE(("Processing: %s" % sinex_name) , level=1)

    ll=sinex_name.split('/')
    sinex_basename=ll[-1]
    
    # if type_glx generates ssc and creates a daily apr file for running glred later
    
    if type_sol == 'glx' :
        glx=sinex_name
        snx=sinex_basename.split('.')[-2]+'.snx'
        dir_snx=path_sinex
        dir_ssc=path_ssc
        ssc=sinex_basename.split('.')[-2]+'.ssc'
        SSinex.glx2snx(glx,snx,dir_snx=dir_snx)
        SSinex.snx2ssc(dir_snx+'/'+snx,ssc,dir_ssc)

        full_sinex_name = dir_snx+'/'+snx
        sinex_name=dir_ssc+'/'+ssc


    ###########################################################################
    # READS FREE SINEX
    ###########################################################################
    
    free_sinex=SSinex.SSinex(name=sinex_name,estimates={})
    free_sinex.read_section_estimate(discontinuity=discontinuities,rename=conf.rename)

    print('----------------------------------------------------------------------------')
    MESSAGE("n_site in free solution: %5d" % len(list(free_sinex.estimates.keys())))

    if args.ref_pck:
        H_common_Gpoint_ref = get_H_common( ref_pck, free_sinex.estimates, prefit=3.5)
        VERBOSE(" Found %03d sites common with %s"  % ( len(H_common_Gpoint_ref), args.ref_pck) )
    else:
        H_common_Gpoint_ref = ref_sinex.common(free_sinex.estimates, prefit=conf.ref_prefit_value)
        VERBOSE("Found %d sites common with %s" % (len(H_common_Gpoint_ref),ref_sinex.name))

    # Check whether Helmert calculation is possible
        
    if len(H_common_Gpoint_ref) < 3:
        WARNING("Number of sites shared with reference solution is less than 3")
        WARNING("Ignoring this date and sinex: %s " % free_sinex.basename)
        T = None
        continue
    
    ###########################################################################
    # HELMERT CALCULATION
    ###########################################################################

    # test if there are excluded components

    if (free_sinex.basename in list(conf.ref_reject.keys())):
        VERBOSE("There are rejected components for %s" % free_sinex.basename)
        lexclude=conf.ref_reject[free_sinex.basename]
    else:
        lexclude={}

    H_common_Gpoint_free=free_sinex.subset_from_keys(H_common_Gpoint_ref)

    # HELMERT

    try:
    
        (T,Residuals,H_stat)=\
            Helmert.estimate_helmert(H_common_Gpoint_ref,H_common_Gpoint_free,vcv_xyz=[],psd=psd,
                            tran=True,rot=True,scale=True,lexclude=lexclude,
                            verbose=verbose,threshold=args.dikin_outlier,method=args.method)
    
    except Exception as error:
        WARNING("Problem in Helmert transformation. skipping %s from solution" % ( free_sinex.basename ) )
        # go to next sinex
        continue
        
    # Rejected components in Helmert estimation have a 1 in column (5,6,7) (east,north,up)
    lindex = np.argwhere(   \
                            ( np.array (list(map(float,Residuals[:,5]))) + \
                             np.array (list(map(float,Residuals[:,6])))  \
                             #np.array (list(map(float,Residuals[:,7])))  \
                            ) == 0 )

    # Get a list of code for common sites used in Helmert transformation ; They will be used for projection of covariance 
    # lcode_proj = np.array(  list( H_common_Gpoint_free.keys()) )[:,0].tolist()
    lcode_proj = Residuals[lindex,0].flatten().tolist()

    MESSAGE("Helmert: n_site_ref: %5d n_obs: %5d n_rejected: %5d   (%5.1lf %%) " % \
            (H_stat['n_site'],H_stat['n_obs_raw'],H_stat['n_obs_rejected'],H_stat['percent_rejected']))

    MESSAGE("Postfit L1 norm median: %10.2lf %10.2lf %10.2lf " % tuple(H_stat['median_L1']))
    MESSAGE("Postfit L2 norm wrms  : %10.2lf %10.2lf %10.2lf " % tuple(H_stat['wrms']))
    MESSAGE("--------------------------------------------------------------------------")



    ###########################################################################
    # UPDATING REFERENCE SITE RESIDUAL TIME SERIES
    ###########################################################################

    # get the date ; all sites will have the same date
    
    ( code, soln ) = list(free_sinex.estimates.keys())[0]
    M=free_sinex.estimates[code,soln]
    
    date = M.epoch

    SORTED_RESIDUALS=Residuals[np.argsort(Residuals[:, 0])]
    
    for i in np.arange(SORTED_RESIDUALS.shape[0]):
        
        code,soln=SORTED_RESIDUALS[i,:2]
        [Re,Rn,Ru] = list(map(float,SORTED_RESIDUALS[i,2:5])) 

        Re = Re * 1.E-3
        Rn = Rn * 1.E-3
        Ru = Ru * 1.E-3


        if l_res_ts.has_ts(code):
            l_res_ts.__dict__[code].add_obs(date,[Rn,Re,Ru,1.E-3,1.E-3,1.E-3,0.,0.,0.],in_place=True,check=False,verbose=False)
        else:
            new_ts = Gts( code = code )
            new_ts.add_obs(date,[Rn,Re,Ru,1.E-3,1.E-3,1.E-3,0.,0.,0.],in_place=True,check=False,verbose=False)
            l_res_ts.append(new_ts)

    H_res_helmert[sinex_basename] =[date]+list(H_stat['wrms'])
    
    ###########################################################################
    # WRITING HELMERT SUMMARY AND STORING TSG VALUE
    ###########################################################################

    VERBOSE("writing Helmert summary")

    if isinstance(T,np.ndarray):
        pyacs.sol.log.helmert_residuals(Residuals, H_stat,free_sinex.name,sum_helmert_file,verbose=verbose)

    TSG[index_sinex,1:]=T.flatten()
    TSG[index_sinex,0]=free_sinex.epoch
    
    l_ref_code += Residuals[:,0].tolist()

    ###########################################################################
    # COORDINATES UNCERTAINTIES
    # change by JMN 01/10/2020
    # if a site has been renamed from a rename command in a conf file
    # then use the old name to get the COV since Paul Rebischung's routine
    # has read the original code

    H_rename = {}

    for code, soln in list(free_sinex.estimates.keys()):
        H_rename[code] = code

    if conf.rename is not None:

        # Case for a CODE rename applying for all SINEX files
        if 'all' in conf.rename:

            for (ccode, new_code) in conf.rename['all']:
                H_rename[ccode] = new_code

        # Case for a CODE rename applying for the current SINEX

        if free_sinex.name in list(conf.rename.keys()):

            for (ccode, new_code) in conf.rename[free_sinex.name]:
                H_rename[ccode] = new_code

    ###########################################################################

    if ( H_COV is None ) and ( args.uncertainty ) and ( type_sol in ['snx','glx'] ):
        VERBOSE("calculating posterior coordinates uncertainties for %s " % ( sinex_name ) )

        from pyacs.sinex.sinex import sinex as pr_sinex
        
        if type_sol == 'glx':
            wsnx_name = full_sinex_name
        else:
            wsnx_name = sinex_name
            
        
        wsnx = pr_sinex.read( wsnx_name )
        lcode = []
        lpt   = []
        lsoln = []
        
        # fills lcode, lpt and lsoln
        for sta in wsnx.sta:
            lcode.append( H_rename[sta.code] )
            lpt.append(sta.pt)
            lsoln.append(sta.soln)

        # unconstrain and make proj
        wsnx.unconstrain()
        wsnx.neqinv(keep_const_on=['UT'], min_const_on='RST', min_const_sig=1e-5 , code = lcode_proj )
        
        H_COV= dict(zip(lcode , wsnx.get_cov_sta(lcode) ) )


    
    ###########################################################################
    # WRITES TRANSFORMED COORDINATES
    ###########################################################################

    VERBOSE("writing transformed coordinates")
        
    snx_transformed=free_sinex.apply_helmert(T, verbose=False)

    ###########################################################################
    # UPDATE TIME SERIES
    ###########################################################################

    if verbose:
        print("-- updating time series from solution: %s" % (sinex_name) )

    dec_date = snx_transformed.epoch
    for code, soln in sorted( snx_transformed.estimates.keys() ):
        M = snx_transformed.estimates[code, soln]

        # new gts ?
        
        if not sts.has_ts( M.code ):
            sts.append( Gts( code = M.code) )
        
        # xyz obs to be added to the time series
        
        XYZSXSYSZCXYCXZCYZ = np.array([ M.X, M.Y, M.Z, 0.001, 0.001, 0.001, 0., 0., 0. ])

        #  uncertainties            
        if H_COV is not None:

            COV = H_COV[code]
            (corr, sigma) = pyacs.lib.glinalg.cov_to_corr(COV)
            (Sx, Sy, Sz) = sigma.flatten()
            Rxy = corr[0, 1]
            Rxz = corr[0, 2]
            Ryz = corr[1, 2]

            XYZSXSYSZCXYCXZCYZ = np.array([ M.X, M.Y, M.Z, Sx, Sy, Sz, Rxy, Rxz, Ryz ])

        # update gts
        sts.__dict__[ M.code ].add_obs_xyz(dec_date,XYZSXSYSZCXYCXZCYZ,in_place=True,check=False, neu=False , verbose=False)
            
    ###########################################################################
    # CLEAN SINEX
    ###########################################################################

    if clean and type_sol == 'glx':
        os.remove(full_sinex_name )

    ###################################################################
    # GLRED SOLUTION IF GLX
    ###################################################################

    if glred_opt and ( type_sol == 'glx' ) :
        VERBOSE("Creating a daily apr file for glred")
            
        apr_4_glred=open(wdir + '/pyacs.apr','w+')
        
        for M in list(snx_transformed.estimates.values()):
            apr_4_glred.write(' %s_GPS %16.4lf  %16.4lf  %16.4lf    0.0000 0.0000 0.0000 %8.4lf 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 from pyacs %s\n' % \
                              (M.code.upper(),M.X,M.Y,M.Z,M.epoch,wdir))

        apr_4_glred.close()

        MESSAGE("Running Glred for %s" % glx)
        SSinex.glred_1_day(glx,path_prt,path_log,wdir + '/pyacs.apr',eq_rename=args.eq_rename)
        os.remove(wdir + '/pyacs.apr')

    index_sinex=index_sinex+1

    MESSAGE(("Processing time since start: %.1lf seconds =  %.1lf minutes " % (time.time()-start,(time.time()-start)/60.)),level=1)



###################################################################
# END LOOP FOR EACH SINEX FILES
###################################################################
MESSAGE("END OF PROCESSING - LOOP ON SINEX FILES",level=1)

###################################################################
# UPDATE STS
###################################################################

if replace:
    MESSAGE("reading original time series")
    try:
        ots = Sgts(ts_dir=path_pck, read=True, verbose=verbose)
    except:
        ots = Sgts(ts_dir=path_pos, read=True, verbose=verbose)

    MESSAGE("%d time series read" % ots.n() )
    MESSAGE("merging original and new time series")

    i = 1
    for code in sorted( sts.lcode() ):
        MESSAGE("merging time series for %s %04d / %04d " % (code, i, len( sts.lcode() )))
        # insert the created ts
        if ots.has_ts(code):
            ots.__dict__[code] = ots.__dict__[code].insert_ts(sts.__dict__[code], rounding='hour',data='xyz',overlap=True)
        # calculates .data from .data_xyz
        ots.__dict__[code].xyz2neu(corr=True)
        i=i+1

    sts = ots.copy()

else:
    MESSAGE("reordering date in times series")
    sts.gts('reorder')

###################################################################
# WRITE TIME SERIES POS FORMAT
###################################################################
if not pck_only:
    MESSAGE("writing time series as pos files")

    #TS=Sgts(path_txyz,verbose=verbose)

    i = 1
    for code in sorted( sts.lcode() ):
        ts = sts.__dict__[code]
        MESSAGE("writing time series for %s %04d / %04d " % (ts.code, i,len( sts.lcode() )))
        ts.xyz2neu(corr=True)
        ts.reorder()
        ts.correct_duplicated_dates(action='correct',tol= .1, in_place=True,verbose=False)
        ts.write_pos(path_pos)
        i = i + 1

###################################################################
# WRITE TIME SERIES PYACS PCK FORMAT
###################################################################

MESSAGE("writing time series as pck file: %s" % (path_pck+'/'+args.expt+'.pck') )

sts.write_pck( path_pck+'/'+args.expt+'.pck', verbose=verbose)

###################################################################
# WRITE TSG TIME SERIES
###################################################################

if not pck_only:
    VERBOSE("writing Helmert parameters values in %s " % path_tsg+'/tsg.dat')

    pyacs.sol.log.write_tsg(lsinex, TSG, path_tsg+'/tsg.dat')

###################################################################
# WRITE RESIDUAL TIME SERIES FOR REFERENCE SITES
###################################################################
VERBOSE("writing residual time series for reference sites in: %s" % path_res_pos)

l_res_ts.gts('write_cats',path_res_pos)

###################################################################
# WRITE STATS FOR RESIDUAL TIME SERIES FOR REFERENCE SITES
###################################################################
VERBOSE("writing statistics for residual time series of reference sites in: %s" %  ( wdir + '/stat/res_ref_stat.dat') )

# get np_code
np_code = np.array( sorted(l_res_ts.lcode()) , dtype=str )

# stat array
stat_array = np.zeros(( np_code.shape[0] , 14 ))
for i in np.arange( np_code.shape[0] ):
    code = np_code[i]
    data = l_res_ts.__dict__[code].data[:,1:] *1.E3
    stat_array[i,0] = data.shape[0]
    stat_array[i,1] = data.shape[0] / len(lsinex) * 100.
    stat_array[i,2:5] = np.mean( data[:,:3] , axis=0 )
    stat_array[i,5:8] = np.median( data[:,:3] , axis=0 )
    stat_array[i,8:11] = np.sqrt(np.mean( data[:,:3]**2 , axis=0 ))
    stat_array[i,11:14] = np.mean(np.fabs(np.diff(data[:,:3] , axis=0 )),axis=0)

header = "Statistics for reference sites\n"
header += "# n: number of occurrence in Helmert parameter estimates\n"
header += "# %: percentage of occurrence in Helmert parameter estimates\n"
header += "# BN,BE,BU: mean bias for NEU in mm\n"
header += "# MBN,MBE,MBU: median bias for NEU in mm\n"
header += "# RN,RE,RU: rms for NEU in mm\n"
header += "# DN,DE,DU: mean daily scatter for NEU in mm\n"
header += "# n    %    BN    BE   BU   MN    ME   MU   RN    RE   RU   DN    DE   DU   \n"
header += "#---------------------------------------------------------------------------\n"

import pyacs.lib.utils
fmt = "%04d %6.1lf %9.2lf %9.2lf %9.2lf %9.2lf %9.2lf %9.2lf %9.2lf %9.2lf %9.2lf %9.2lf %9.2lf %9.2lf %s"
pyacs.lib.utils.save_np_array_with_string(stat_array,np_code, fmt, wdir + '/stat/res_ref_stat.dat' , header)

###################################################################
# COMMAND LINE
###################################################################

cmd_line_file=open(wdir + '/stat/pyacs_cmd_line.dat','w')
VERBOSE("writing command line in: %s" % ( wdir + '/stat/pyacs_cmd_line.dat') )

cmd_line_file.write("%s" % cmd_line)
cmd_line_file.close()

###################################################################
# HELMERT WRMS TIME SERIES
###################################################################

htswrms_file=open(wdir + '/stat/helmert_wrms_sum.dat','w')
VERBOSE("writing helmert wrms time series in : %s" % (wdir + '/stat/helmert_wrms_sum.dat') )

for key,value in H_res_helmert.items():
    htswrms_file.write("%-40s %9.5lf %8.2lf %8.2lf %8.2lf\n" % tuple([key]+list(value)))

htswrms_file.write("#--------------------------------------------------------------------------------\n")
htswrms_file.write("#%48s   %8.2lf %8.2lf %8.2lf\n" % tuple( ["mean (East, North, Up, in mm)"] +list(np.mean(np.array(list(H_res_helmert.values())),axis=0)[1:]) ))

htswrms_file.close()


###################################################################
# APR
###################################################################
if not pck_only:
    apr_file= path_stat + '/pyacs.apr'
    MESSAGE("writing apr file in %s" % apr_file)

    i = 1
    for code in sorted( sts.lcode() ):
        VERBOSE("%s %04d / %04d" % (code,i,len( sts.lcode() )))
        ts = sts.__dict__[code]
        ts.save_apr(apr_file, epoch=None, verbose=False, excluded_periods=None)
        i = i + 1
    
###################################################################
# INFO FILES
###################################################################
if not pck_only:
    if sts.lcode() == []:
        WARNING("info file empty. no site processed")

    else:
        MESSAGE("writing information files in %s" % path_stat)

        sts.stat_site(save_dir=path_stat)

        ###################################################################
        # VELOCITY MAPS GMT PSVELO AND SHAPEFILES
        ###################################################################

        import pyacs.lib.shapefile

        src = path_stat+'/pyacs_vel.dat'
        dest= path_maps+'/pyacs_vel.dat'
        shutil.move(src, dest)

        gmt = dest
        shp = path_maps+'/pyacs_vel.shp'
        pyacs.lib.shapefile.psvelo_to_shapefile(gmt, shp, verbose=verbose)

        src = path_stat+'/pyacs_vel_up.dat'
        dest= path_maps+'/pyacs_vel_up.dat'
        shutil.move(src, dest)

        shp = path_maps+'/pyacs_vel_up.shp'
        pyacs.lib.shapefile.psvelo_to_shapefile(gmt, shp, verbose=verbose)


###############################################################################
# CREATES GLRED POS FILE IF TYPE_GLX
###############################################################################
if not pck_only:
    if type_sol == 'glx' and glred_opt :

        MESSAGE("Creating PBO format pos files")
        import subprocess
        os.chdir(path_prt)
        cmd='sh_plot_pos -f prt.org -p '
        MESSAGE("Running %s" % cmd)
        subprocess.getstatusoutput(cmd)
        lpos=glob.glob('*.pos')
        for pos in lpos:
            shutil.move(pos, '../glred_pos')

###############################################################################
# REMOVES USELESS DIRECTORIES
###############################################################################

if not glred_opt :
    VERBOSE("Removing glred directories of %s " % args.expt)
    try:
        shutil.rmtree(path_lgred_pos,ignore_errors=True)
        shutil.rmtree(path_prt,ignore_errors=True)
        shutil.rmtree(path_log,ignore_errors=True)
        
    except:
        VERBOSE("Problem Removing glred directories of %s" % args.expt)
        
        pass
    
if clean:
    VERBOSE("Removing glred_sinex directory of %s" % args.expt)
    try:
        shutil.rmtree(path_sinex)
        
    except:
        pass


###################################################################
# WRITE RUNTIME
###################################################################
MESSAGE(("-- Processing time since start: %.1lf seconds =  %.1lf minutes " % (time.time()-start,(time.time()-start)/60.)),level=1)
