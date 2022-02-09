#-------------------------------------------------------------------------------
# Module   : sinex
# Purpose  : Read, write and manipulate SINEX files
# Author   : P. Rebischung
# Created  : 20-May-2011
#
# Changes  :
#
# Routines : - __init__            : Initialize a sinex object
#            - read                : Read a SINEX file into a sinex object
#            - write               : Write a sinex object into a SINEX file
#            - check_staid         : Check PT codes and DOMES numbers
#            - check_solns         : Check solns in an "instantaneous" solution
#            - check_erp           : Set ERP reference epochs to exactly noon
#            - get_xyz             : Get cartesian coordinates of specified stations
#            - get_plh             : Get geographical coordinates of specified stations
#            - get_lonlat          : Get longitudes and latitudes of specified stations
#            - get_sigenh          : Get E,N,H formal errors of specified stations
#            - get_helmert_matrix  : Get design matrix of Helmert parameters
#            - map                 : Draw station map
#            - propagate           : Propagate station positions to specified date
#            - remove_sta          : Remove specified stations from a solution
#            - keep_sta            : Keep only specified stations in a solution
#            - remove_sta_wo_domes : Remove stations without DOMES numbers from a solution
#            - keep_relevant_solns : Extract solns that a relevant to specified epoch from a solution
#            - remove_params       : Remove certain types of parameters from a solution
#            - remove_undef_params : Remove "undefined" parameters from a solution
#            - unconstrain         : Recover unconstrained normal equation from a solution
#            - prior2ref           : Change a priori values to reference values in a normal equation
#            - fix_params          : Fix certain types of parameters in a normal equation
#            - reduce_params       : Reduce certain types of parameters from a normal equation
#            - reduce_sta          : Reduce specified stations from a normal equation
#            - reduce_doublons     : Reduce stations appearing twice in a normal equation
#            - setup_geocenter     : Add geocenter coordinates into a normal equation
#            - setup_scale         : Add terrestrial scale into a normal equation
#            - add_min_const       : Add minimal constraints to normal matrix of constraints
#            - set_constraints     : Define normal matrix of constraints
#            - neqinv              : Invert normal equation
#            - rescale             : Multiply solution covariance matrix by given factor
#            - get_common_points   : Get list of common points between two solutions
#            - helmert_wrt         : Helmert comparison between two solutions
#            - trim_metadata       : Remove metadata that are not relevant for specified period
#            - update_stalist      : Given an AC solution, update list of stations in the IGS combined solution
#            - get_core            : Get list of core stations in a solution
#            - extract_core        : Extract core stations from a solution
#            - calib_lod           : Calibrate LOD estimates wrt reference series (Bulletin A)
#            - apriori_sf          : Compute a priori scale factor of a solution
#            - add_metadata        : Add metadata blocks into a sinex object
#            - get_residuals       : Compare solution to a reference solution and get necessary information for res file
#            - check_sat_pco       : Check satellite antenna phase center offsets
#            - correct_opoleload   : Correct ocean pole tide loading displacements in a normal equation
#            - write_solndomes     : Write special soln file with DOMES numbers (needed for IGS combination database)
#            - get_nonpubsolns     : Get list of solns to remove from official IGS cumulative solution
#            - add_psd             : Add post-seismic deformation models to a solution
#            - remove_psd          : Remove post-seismic deformation models from a solution
#            - add_per             : Add periodic terms to a solution
#            - plate_poles         : Estimate tectonic plate rotation poles from site velocities
#            - insert_disc         : Insert new "discontinuity" for a given station, by duplicating appropriate soln
#            - apply_offsets       : Apply position offsets to specified solns
#            - write_solndomes     : Write special soln file with DOMES numbers
#            - loadest_wrt         : Estimate surface load coefficients from comparison to a reference solution
#-------------------------------------------------------------------------------



# LIBRARIES
#-------------------------------------------------------------------------------
import copy
from math import *
import numpy
from scipy import sparse
from scipy import special

from .constants import *
from .utils import *
from .snxutils import *
from .mathutils import *
from .date import date



# CONSTANTS
#-------------------------------------------------------------------------------

# List of blocks of the SINEX format (unsupported blocks are commented)
blocks = ['FILE/REFERENCE',                   #  0
          'FILE/COMMENT',                     #  1
          'INPUT/HISTORY',                    #  2
          'INPUT/FILES',                      #  3
          'INPUT/ACKNOWLEDGEMENTS',           #  4
          'SITE/ID',                          #  5
          'SITE/RECEIVER',                    #  6
          'SITE/ANTENNA',                     #  7
          'SITE/GPS_PHASE_CENTER',            #  8
          'SITE/ECCENTRICITY',                #  9
          'SOLUTION/EPOCHS',                  # 10
          'SOLUTION/ESTIMATE',                # 11
          'SOLUTION/APRIORI',                 # 12
          'SOLUTION/MATRIX_ESTIMATE',         # 13
          'SOLUTION/MATRIX_APRIORI',          # 14
          'SOLUTION/NORMAL_EQUATION_VECTOR',  # 15
          'SOLUTION/NORMAL_EQUATION_MATRIX',  # 16
          'SOLUTION/STATISTICS']              # 17
          #'NUTATION/DATA',
          #'PRECESSION/DATA',
          #'SOURCE/ID',
          #'SITE/DATA',
          #'SITE/GAL_PHASE_CENTER',
          #'SATELLITE/ID',
          #'SATELLITE/PHASE_CENTER',
          #'BIAS/EPOCHS']

# Default DOMES number
default_domes = '     M   '



# SINEX CLASS
#-------------------------------------------------------------------------------
class sinex:

#-------------------------------------------------------------------------------
# Routine : __init__
# Purpose : Initialize a sinex object
# Author  : P. Rebischung
# Created : 20-May-2011
#
# Changes :
#
# Input   :
# Output  :
#-------------------------------------------------------------------------------
    def __init__(self):

        self.file     = None
        self.version  = None
        self.agency   = None
        self.epoch    = None
        self.start    = None
        self.end      = None
        self.tech     = None
        self.npar     = None
        self.const    = None
        self.content  = None
        self.ref      = None
        self.comment  = None
        self.input    = None
        self.acks     = None
        self.sta      = None
        self.antenna  = None
        self.param    = None
        self.x        = None
        self.sig      = None
        self.prior    = None
        self.x0       = None
        self.sig0     = None
        self.Q        = None
        self.Ntot     = None
        self.Qc       = None
        self.Nc       = None
        self.k        = None
        self.N        = None
        self.nobs     = None
        self.nunk     = None
        self.sampling = None
        self.vPv      = None
        self.sigphase = None
        self.sigcode  = None
        self.dof      = None
        self.vf       = None
        self.lPl      = None

#-------------------------------------------------------------------------------
# Routine : read
# Purpose : Read a SINEX file into a sinex object
# Author  : P. Rebischung
# Created : 20-May-2011
#
# Changes : PR, 04-Sep-2014 : Read SOLUTION/STATISTICS block
#
# Input   : - file : Name of SINEX file to read
#           - dont_read : List of keywords to indicate which blocks should not
#             be read. dont_read can include the following keywords:
#             'matrices' in order not to read matrices
#             'comments' in order not to read all "comment" blocks
#             'apriori'  in order not to read a priori information
#             'metadata' in order not to read receivers, antennas...
# Output  : - snx : sinex object
#-------------------------------------------------------------------------------
    @classmethod
    def read(cls, file, dont_read = []):

        # Initializations
        snx = sinex()
        snx.file = file
        adress = [0]*len(blocks)
        read   = [1]*len(blocks)
        type_matest = ''
        type_matapr = ''

        # In function of argument dont_read, mark blocks which should not be read.
        if ('matrices' in dont_read):
            read[13] = 0
            read[14] = 0
            read[16] = 0
        if ('comments' in dont_read):
            read[0:5] = [0]*5
        if ('apriori' in dont_read):
            read[12] = 0
            read[14] = 0
        if ('metadata' in dont_read):
            read[6:10] = [0]*4

        # Open input SINEX file
        f = open(file, 'r')

        # Read 1st line
        line = f.readline()
        snx.version = line[6:10]
        snx.agency = line[11:14]

# Don't know why does work any more in Python3.6        
#        snx.epoch = re.sub(' ', '0', line[15:27])
#        snx.start = re.sub(' ', '0', line[32:44])
#        snx.end = re.sub(' ', '0', line[45:57])
        #print(line[15:27])
        #print(type(line[15:27]))
        snx.epoch = line[15:27].replace(' ','0')
        snx.start = line[32:44].replace(' ','0')
        snx.end   = line[45:57].replace(' ','0')
        
        snx.tech = line[58:59]
        snx.const = line[66:67]
        snx.content = line[68:].strip()

        # Read rest of the file to get adress of the blocks
        # Also get types of the matrices (COVA/CORR/INFO)
        while (line):
            if (line[0:1] == '+'):
                block = line[1:line.find(' ')].strip()
                if (block in blocks):
                    i = blocks.index(block)
                    adress[i] = f.tell()
                if (line[1:25] == 'SOLUTION/MATRIX_ESTIMATE'):
                    type_matest = line[28:32]
                elif (line[1:24] == 'SOLUTION/MATRIX_APRIORI'):
                    type_matapr = line[27:31]
            line = f.readline()

        # Read FILE/REFERENCE block -> snx.ref
        if (adress[0] and read[0]):
            snx.ref = record()
            snx.ref.description = None
            snx.ref.output      = None
            snx.ref.contact     = None
            snx.ref.software    = None
            snx.ref.hardware    = None
            snx.ref.input       = None
            f.seek(adress[0])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    if (line[1:12] == 'DESCRIPTION'):
                        snx.ref.description = line[20:].strip()
                    elif (line[1:7] == 'OUTPUT'):
                        snx.ref.output = line[20:].strip()
                    elif (line[1:8] == 'CONTACT'):
                        snx.ref.contact = line[20:].strip()
                    elif (line[1:9] == 'SOFTWARE'):
                        snx.ref.software = line[20:].strip()
                    elif (line[1:9] == 'HARDWARE'):
                        snx.ref.hardware = line[20:].strip()
                    elif (line[1:6] == 'INPUT'):
                        snx.ref.input = line[20:].strip()
                line = f.readline()

        # Read FILE/COMMENT block -> snx.comment
        if (adress[1] and read[1]):
            snx.comment = []
            f.seek(adress[1])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    snx.comment.append(line[1:].strip())
                line = f.readline()

        # Read INPUT/HISTORY block -> snx.input
        if (adress[2] and read[2]):
            snx.input = []
            f.seek(adress[2])
            line = f.readline()

            while (line[0:1] != '-'):
                if ((line[0:1] != '*') and (line[1:2] == '+')):
                    r = record()
                    r.version = line[6:10]
                    r.agency  = line[11:14]
                    r.epoch   = line[15:27]
                    r.start   = line[32:44]
                    r.end     = line[45:57]
                    r.tech    = line[58:59]
                    r.npar    = int(line[60:65])
                    r.const   = line[66:67]
                    r.content = line[68:].strip()
                    snx.input.append(r)
                line = f.readline()

        # Read INPUT/FILES block -> complement snx.input
        if (adress[3] and read[3]):
            f.seek(adress[3])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    agency = line[1:4]
                    epoch  = line[5:17]
                    if (agency+epoch in [inp.agency+inp.epoch for inp in snx.input]):
                        i = [inp.agency+inp.epoch for inp in snx.input].index(agency+epoch)
                        snx.input[i].file = line[18:47]
                        snx.input[i].description = line[48:].strip()
                line = f.readline()

        # Read INPUT/ACKNOWLEDGEMENTS block -> snx.acks
        if (adress[4] and read[4]):
            snx.acks = []
            f.seek(adress[4])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    r = record()
                    r.agency = line[1:4]
                    r.description = line[5:].strip()
                    snx.acks.append(r)
                line = f.readline()

        # Read SITE/ID block -> snx.sta
        if (adress[5] and read[5]):
            snx.sta = []
            f.seek(adress[5])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    r = record()
                    r.code         = line[1:5].upper()
                    r.pt           = line[6:8]
                    r.domes        = line[9:18]
                    r.tech         = line[19:20]
                    r.description  = line[21:43]
                    r.lon          = line[44:55]
                    r.lat          = line[56:67]
                    r.h            = line[68:75]
                    r.receiver     = []
                    r.antenna      = []
                    r.eccentricity = []
                    r.soln         = []
                    snx.sta.append(r)
                line = f.readline()

        # Read SITE/RECEIVER block -> snx.sta[*].receiver
        if (adress[6] and read[6]):
            f.seek(adress[6])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    code = line[1:5].upper()
                    pt   = line[6:8]

                    # PATCH: Add possibly missing station into snx.sta
                    if (not(code+pt in [s.code+s.pt for s in snx.sta])):
                        r = record()
                        r.code = code
                        r.pt = pt
                        r.domes = default_domes
                        r.tech = 'P'
                        r.description = 22*' '
                        r.lon = 11*' '
                        r.lat = 11*' '
                        r.h = 7*' '
                        r.receiver = []
                        r.antenna = []
                        r.eccentricity = []
                        r.soln = []
                        snx.sta.append(r)

                    i = [s.code+s.pt for s in snx.sta].index(code+pt)
                    r = record()
                    r.start    = line[16:28]
                    r.end      = line[29:41]
                    r.type     = line[42:62].strip()
                    while (len(r.type) < 20):
                        r.type = r.type + ' '
                    r.serie    = line[63:68].strip()
                    if (r.serie == ''):
                        r.serie = '-----'
                    while (len(r.serie) < 5):
                        r.serie = r.serie + ' '
                    r.firmware = line[69:].strip()
                    snx.sta[i].receiver.append(r)
                line = f.readline()

        # Read SITE/ANTENNA block -> snx.sta[*].antenna
        if (adress[7] and read[7]):
            f.seek(adress[7])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    code = line[1:5].upper()
                    pt   = line[6:8]
                    i = [s.code+s.pt for s in snx.sta].index(code+pt)
                    r = record()
                    r.start    = line[16:28]
                    r.end      = line[29:41]
                    r.type     = line[42:62]
                    if (r.type[16:] == '    '):
                        r.type = r.type[0:16]+'NONE'
                    r.serie    = line[63:68].strip()
                    snx.sta[i].antenna.append(r)
                line = f.readline()

        # Read SITE/GPS_PHASE_CENTER block -> snx.antenna
        if (adress[8] and read[8]):
            snx.antenna = []
            f.seek(adress[8])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    r = record()
                    r.type = line[1:21]
                    r.serie = line[22:27]
                    r.pco = [0]*6
                    r.pco[0] = line[28:34]
                    r.pco[1] = line[35:41]
                    r.pco[2] = line[42:48]
                    r.pco[3] = line[49:55]
                    r.pco[4] = line[56:62]
                    r.pco[5] = line[63:69]
                    r.model  = line[70:].strip()
                    snx.antenna.append(r)
                line = f.readline()

        # Read SITE/ECCENTRICITY block -> snx.sta[*].eccentricity
        if (adress[9] and read[9]):
            f.seek(adress[9])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    code = line[1:5].upper()
                    pt   = line[6:8]
                    i = [s.code+s.pt for s in snx.sta].index(code+pt)
                    r = record()
                    r.start    = line[16:28]
                    r.end      = line[29:41]
                    r.system   = line[42:45]
                    r.dx       = [0]*3
                    r.dx[0]    = line[46:54]
                    r.dx[1]    = line[55:63]
                    r.dx[2]    = line[64:72]
                    snx.sta[i].eccentricity.append(r)
                line = f.readline()

        # Read SOLUTION/EPOCHS block -> snx.sta[*].soln
        if (adress[10] and read[10]):
            f.seek(adress[10])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    code = line[1:5].upper()
                    pt   = line[6:8]
                    if (code+pt in [s.code+s.pt for s in snx.sta]):
                        i = [s.code+s.pt for s in snx.sta].index(code+pt)
                        r = record()
                        r.soln = line[9:13]
#                        r.datastart = re.sub(' ', '0', line[16:28])
                        r.datastart = line[16:28].replace(' ','0')
                        if (r.datastart == '00:000:00000'):
                            r.datastart = snx.start
#                        r.dataend = re.sub(' ', '0', line[29:41])
                        r.dataend = line[29:41].replace(' ','0')
                        if (r.dataend == '00:000:00000'):
                            r.dataend = snx.end
#                        r.datamean = re.sub(' ', '0', line[42:54])
                        r.datamean = line[42:54].replace(' ','0')
                        if (r.datamean == '00:000:00000'):
                            mjd1 = date.from_snxepoch(r.datastart).mjd
                            mjd2 = date.from_snxepoch(r.dataend).mjd
                            r.datamean = date.from_mjd((mjd1 + mjd2) / 2).snxepoch()
                        snx.sta[i].soln.append(r)
                line = f.readline()

        # Read SOLUTION/ESTIMATE block -> snx.param, snx.x and snx.sig
        if (adress[11] and read[11]):
            snx.param = []
            snx.x = []
            snx.sig = []
            f.seek(adress[11])
            line = f.readline()

            i = -1
            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    i = i+1
                    r = record()
                    r.type  = line[7:13]
                    r.code  = line[14:18].upper()
                    r.pt    = line[19:21]
                    r.soln  = line[22:26]
#                    r.tref  = re.sub(' ', '0', line[27:39])
                    r.tref  = line[27:39].replace(' ', '0')
                    r.unit  = line[40:44]
                    r.const = line[45:46]
                    snx.param.append(r)
                    snx.x.append(float(line[47:68]))
                    snx.sig.append(float(line[69:80]))
                line = f.readline()

            snx.x = numpy.array(snx.x)
            snx.sig = numpy.array(snx.sig)
            snx.npar = len(snx.param)

        # Read SOLUTION/APRIORI block -> snx.prior, snx.x0 and snx.sig0
        if (adress[12] and read[12]):
            snx.prior = []
            snx.x0 = []
            snx.sig0 = []
            f.seek(adress[12])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    r = record()
                    r.type  = line[7:13]
                    r.code  = line[14:18].upper()
                    r.pt    = line[19:21]
                    r.soln  = line[22:26]
                    r.tref  = line[27:39].replace(' ', '0')
                    r.unit  = line[40:44]
                    r.const = line[45:46]
                    snx.prior.append(r)
                    snx.x0.append(float(line[47:68]))
                    snx.sig0.append(float(line[69:80]))
                line = f.readline()

            snx.x0 = numpy.array(snx.x0)
            snx.sig0 = numpy.array(snx.sig0)

        # Read SOLUTION/MATRIX_ESTIMATE block -> snx.Q or snx.Ntot
        if (adress[13] and read[13]):
            f.seek(adress[13])
            line = f.readline()

            # Read matrix
            Q = numpy.zeros((snx.npar, snx.npar))
            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    i = int(line[1:6]) - 1
                    j = int(line[7:12]) - 1
                    Q[i,j] = float(line[13:34])
                    Q[j,i] = Q[i,j]
                    if (line[35:56].strip()):
                        Q[i,j+1] = float(line[35:56])
                        Q[j+1,i] = Q[i,j+1]
                    if (line[57:78].strip()):
                        Q[i,j+2] = float(line[57:78])
                        Q[j+2,i] = Q[i,j+2]
                line = f.readline()

            # If covariance or correlation matrix, store matrix in snx.Q
            if (type_matest in ['CORR', 'COVA']):
                snx.Q = Q

                # If correlation matrix, compute covariances
                if (type_matest == 'CORR'):
                    for i in range(snx.npar):
                        for j in range(i):
                            snx.Q[i,j] = snx.Q[i,j] * snx.Q[i,i] * snx.Q[j,j]
                            snx.Q[j,i] = snx.Q[i,j]
                    for i in range(snx.npar):
                        snx.Q[i,i] = snx.Q[i,i]**2

            # If normal matrix, store matrix in snx.Ntot
            elif (type_matest == 'INFO'):
                snx.Ntot = Q

        # Read SOLUTION/MATRIX_APRIORI block -> snx.Qc or snx.Nc
        if (adress[14] and read[14]):
            f.seek(adress[14])
            line = f.readline()

            # Read matrix
            Q = numpy.zeros((snx.npar, snx.npar))
            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    i = int(line[1:6]) - 1
                    j = int(line[7:12]) - 1
                    Q[i,j] = float(line[13:34])
                    Q[j,i] = Q[i,j]
                    if (line[35:56].strip()):
                        Q[i,j+1] = float(line[35:56])
                        Q[j+1,i] = Q[i,j+1]
                    if (line[57:78].strip()):
                        Q[i,j+2] = float(line[57:78])
                        Q[j+2,i] = Q[i,j+2]
                line = f.readline()

            # If covariance or correlation matrix, store matrix in snx.Qc
            if (type_matapr in ['CORR', 'COVA']):
                snx.Qc = Q

                # If correlation matrix, compute covariances
                if (type_matapr == 'CORR'):
                    for i in range(snx.npar):
                        for j in range(i):
                            snx.Qc[i,j] = snx.Qc[i,j] * snx.Qc[i,i] * snx.Qc[j,j]
                            snx.Qc[j,i] = snx.Qc[i,j]
                    for i in range(snx.npar):
                        snx.Qc[i,i] = snx.Qc[i,i]**2

            # If normal matrix, store matrix in snx.Nc
            elif (type_matapr == 'INFO'):
                snx.Nc = Q

        # Read SOLUTION/NORMAL_EQUATION_VECTOR block -> snx.k
        if (adress[15] and read[15]):
            snx.k = numpy.zeros(snx.npar)
            f.seek(adress[15])
            line = f.readline()

            i = -1
            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    i = i+1
                    snx.k[i] = float(line[47:68])
                line = f.readline()

        # Read SOLUTION/NORMAL_EQUATION_MATRIX block -> snx.N
        if (adress[16] and read[16]):
            f.seek(adress[16])
            line = f.readline()
            snx.N = numpy.zeros((snx.npar, snx.npar))

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    i = int(line[1:6]) - 1
                    j = int(line[7:12]) - 1
                    snx.N[i,j] = float(line[13:34])
                    snx.N[j,i] = snx.N[i,j]
                    if (line[35:56].strip()):
                        snx.N[i,j+1] = float(line[35:56])
                        snx.N[j+1,i] = snx.N[i,j+1]
                    if (line[57:78].strip()):
                        snx.N[i,j+2] = float(line[57:78])
                        snx.N[j+2,i] = snx.N[i,j+2]
                line = f.readline()

        # Read SOLUTION/STATISTICS block
        if (adress[17] and read[17]):
            f.seek(adress[17])
            line = f.readline()

            while (line[0:1] != '-'):
                if (line[0:1] != '*'):
                    if (line[1:31] == 'NUMBER OF OBSERVATIONS        '):
                        snx.nobs = float(line[32:])
                    elif (line[1:31] == 'NUMBER OF UNKNOWNS            '):
                        snx.nunk = float(line[32:])
                    elif (line[1:31] == 'SAMPLING INTERVAL (SECONDS)   '):
                        snx.sampling = float(line[32:])
                    elif (line[1:31] == 'SQUARE SUM OF RESIDUALS (VTPV)'):
                        snx.vPv = float(line[32:])
                    elif (line[1:31] == 'PHASE MEASUREMENTS SIGMA      '):
                        snx.sigphase = float(line[32:])
                    elif (line[1:31] == 'CODE MEASUREMENTS SIGMA       '):
                        snx.sigcode = float(line[32:])
                    elif (line[1:31] == 'NUMBER OF DEGREES OF FREEDOM  '):
                        snx.dof = float(line[32:])
                    elif (line[1:31] == 'VARIANCE FACTOR               '):
                        snx.vf = float(line[32:])
                    elif (line[1:31] == 'WEIGHTED SQUARE SUM OF O-C    '):
                        snx.lPl = float(line[32:])
                line = f.readline()

        # Close input SINEX file
        f.close()

        return snx

#-------------------------------------------------------------------------------
# Routine : write
# Purpose : Write a sinex object into a SINEX file
# Author  : P. Rebischung
# Created : 22-May-2011
#
# Changes :
#
# Input   : - file         : Name of SINEX file to write
#           - dont_write   : List of keywords to indicate which blocks should not
#                            be written. dont_write can include the following keywords:
#                            'matrices' in order not to write matrices
#                            'comments' in order not to write all "comment" blocks
#                            'apriori'  in order not to write a priori information
#                            'metadata' in order not to write receivers, antennas...
#-------------------------------------------------------------------------------
    def write(self, file, dont_write=[], epochs_style='new'):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        # Open output SINEX file
        f = open(file, 'w')

        # Update snx.epoch and write first line
        snx.epoch = date().snxepoch()
        f.write('%=SNX 2.02 {0} {1} {0} '.format(snx.agency, snx.epoch))
        f.write('{0} {1} {2} '.format(snx.start, snx.end, snx.tech))
        f.write('{0:>5} {1} {2}\n'.format(snx.npar, snx.const, snx.content))

        # Write FILE/REFERENCE block
        if ((snx.ref) and (not('comments' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+FILE/REFERENCE\n')
            if (snx.ref.description):
                f.write(' {0:<18} {1}\n'.format('DESCRIPTION', snx.ref.description))
            if (snx.ref.output):
                f.write(' {0:<18} {1}\n'.format('OUTPUT', snx.ref.output))
            if (snx.ref.contact):
                f.write(' {0:<18} {1}\n'.format('CONTACT', snx.ref.contact))
            if (snx.ref.software):
                f.write(' {0:<18} {1}\n'.format('SOFTWARE', snx.ref.software))
            if (snx.ref.hardware):
                f.write(' {0:<18} {1}\n'.format('HARDWARE', snx.ref.hardware))
            if (snx.ref.input):
                f.write(' {0:<18} {1}\n'.format('INPUT', snx.ref.input))
            f.write('-FILE/REFERENCE\n')

        # Write FILE/COMMENT block
        if ((snx.comment) and (not('comments' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+FILE/COMMENT\n')
            for c in snx.comment:
                f.write(' {0}\n'.format(c))
            f.write('-FILE/COMMENT\n')

        # Write INPUT/ACKNOWLEDGEMENTS block
        if ((snx.acks) and (not('comments' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+INPUT/ACKNOWLEDGEMENTS\n')
            f.write('*AGY ______________________________FULL_DESCRIPTION_____________________________\n')
            for a in snx.acks:
                f.write(' {0} {1}\n'.format(a.agency, a.description))
            f.write('-INPUT/ACKNOWLEDGEMENTS\n')

        # Write INPUT/HISTORY block
        if ((snx.input) and (not('comments' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+INPUT/HISTORY\n')
            f.write('*_VERSION_ CRE __CREATION__ OWN _DATA_START_ __DATA_END__ T PARAM S ____TYPE____\n')
            for i in snx.input:
                f.write(' +SNX {0} {1} {2} {1} '.format(i.version, i.agency, i.epoch))
                f.write('{0} {1} {2} '.format(i.start, i.end, i.tech))
                f.write('{0:>5} {1} {2}\n'.format(i.npar, i.const, i.content))
            f.write(' =SNX 2.02 {0} {1} {0} '.format(snx.agency, snx.epoch))
            f.write('{0} {1} {2} '.format(snx.start, snx.end, snx.tech))
            f.write('{0:>5} {1} {2}\n'.format(snx.npar, snx.const, snx.content))
            f.write('-INPUT/HISTORY\n')

        # Write INPUT/FILES block
        if ((snx.input) and (not('comments' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+INPUT/FILES\n')
            f.write('*OWN __CREATION__ ___________FILENAME__________ ___________DESCRIPTION__________\n')
            for i in snx.input:
                f.write(' {0} {1} {2} {3}\n'.format(i.agency, i.epoch, i.file, i.description))
            f.write('-INPUT/FILES\n')

        # Write SITE/ID block
        if (snx.sta):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SITE/ID\n')
            f.write('*CODE PT __DOMES__ T _STATION DESCRIPTION__ _LONGITUDE_ _LATITUDE__ HEIGHT_\n')
            for s in snx.sta:
                f.write(' {0} {1} {2} {3} '.format(s.code, s.pt, s.domes, s.tech))
                f.write('{0} {1} {2} {3}\n'.format(s.description, s.lon, s.lat, s.h))
            f.write('-SITE/ID\n')

        # Write SITE/RECEIVER block
        if ((snx.sta) and (not('metadata' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SITE/RECEIVER\n')
            f.write('*CODE PT SOLN T _DATA START_ __DATA_END__ ___RECEIVER_TYPE____ _S/N_ _FIRMWARE__\n')
            for s in snx.sta:
                for r in s.receiver:
                    f.write(' {0} {1} ---- {2} '.format(s.code, s.pt, s.tech))
                    f.write('{0} {1} {2} {3} {4}\n'.format(r.start, r.end, r.type, r.serie, r.firmware))
            f.write('-SITE/RECEIVER\n')

        # Write SITE/ANTENNA block
        if ((snx.sta) and (not('metadata' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SITE/ANTENNA\n')
            f.write('*CODE PT SOLN T _DATA START_ __DATA_END__ ____ANTENNA_TYPE____ _S/N_\n')
            for s in snx.sta:
                for a in s.antenna:
                    f.write(' {0} {1} ---- {2} '.format(s.code, s.pt, s.tech))
                    f.write('{0} {1} {2} {3}\n'.format(a.start, a.end, a.type, a.serie))
            f.write('-SITE/ANTENNA\n')

        # Write SITE/GPS_PHASE_CENTER block
        if ((snx.antenna) and (not('metadata' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SITE/GPS_PHASE_CENTER\n')
            f.write('*________TYPE________ _S/N_ _L1_U_ _L1_N_ _L1_E_ _L2_U_ _L2_N_ _L2_E_ __MODEL___\n')
            for a in snx.antenna:
                f.write(' {0} {1} '.format(a.type, a.serie))
                for i in range(6):
                    f.write('{0} '.format(a.pco[i]))
                f.write('{0}\n'.format(a.model))
            f.write('-SITE/GPS_PHASE_CENTER\n')

        # Write SITE/ECCENTRICITY block
        if ((snx.sta) and (not('metadata' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SITE/ECCENTRICITY\n')
            f.write('*CODE PT SOLN T _DATA START_ __DATA_END__ REF __DX_U__ __DX_N__ __DX_E__\n')
            for s in snx.sta:
                for e in s.eccentricity:
                    f.write(' {0} {1} ---- {2} '.format(s.code, s.pt, s.tech))
                    f.write('{0} {1} {2} '.format(e.start, e.end, e.system))
                    f.write('{0} {1} {2}\n'.format(e.dx[0], e.dx[1], e.dx[2]))
            f.write('-SITE/ECCENTRICITY\n')

        # Write SOLUTION/EPOCHS block
        if (snx.sta):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/EPOCHS\n')
            f.write('*CODE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_\n')
            for s in snx.sta:
                for i in s.soln:
                    f.write(' {0} {1} {2} {3} '.format(s.code, s.pt, i.soln, s.tech))
                    f.write('{0} {1} {2}\n'.format(i.datastart, i.dataend, i.datamean))
            f.write('-SOLUTION/EPOCHS\n')

        # Write SOLUTION/APRIORI block
        if ((snx.prior) and (not('apriori' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/APRIORI\n')
            f.write('*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ____APRIORI_VALUE____ __STD_DEV__\n')
            for i in range(len(snx.prior)):
                p = snx.prior[i]
                f.write(' {0:5} {1} {2} {3} {4} '.format(i+1, p.type, p.code, p.pt, p.soln))
                f.write('{0} {1} {2} {3:21.14e} {4:11.5e}\n'.format(p.tref, p.unit, p.const, snx.x0[i], snx.sig0[i]))
            f.write('-SOLUTION/APRIORI\n')

        # Write SOLUTION/ESTIMATE block
        if (snx.param):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/ESTIMATE\n')
            f.write('*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__\n')
            for i in range(snx.npar):
                p = snx.param[i]
                f.write(' {0:5} {1} {2} {3} {4} '.format(i+1, p.type, p.code, p.pt, p.soln))
                f.write('{0} {1} {2} {3:21.14e} {4:11.5e}\n'.format(p.tref, p.unit, p.const, snx.x[i], snx.sig[i]))
            f.write('-SOLUTION/ESTIMATE\n')

        # Write SOLUTION/MATRIX_APRIORI block (case of covariance matrix)
        if ((snx.Qc != None) and (not('matrices' in dont_write)) and (not('apriori' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/MATRIX_APRIORI L COVA\n')
            f.write('*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______\n')
            for i in range(len(snx.prior)):
                j = 0
                while (j <= i):
                    if (snx.Qc[i,j] != 0):
                        f.write(' {0:5} {1:5} {2:21.14e}'.format(i+1, j+1, snx.Qc[i,j]))
                        if ((j+2 <= i) and (snx.Qc[i,j+2] != 0)):
                            f.write(' {0:21.14e}'.format(snx.Qc[i,j+1]))
                            f.write(' {0:21.14e}\n'.format(snx.Qc[i,j+2]))
                            j = j+3
                        elif ((j+1 <= i) and (snx.Qc[i,j+1] != 0)):
                            f.write(' {0:21.14e}\n'.format(snx.Qc[i,j+1]))
                            j = j+2
                        else:
                            f.write('\n')
                            j = j+1
                    else:
                        j = j+1
            f.write('-SOLUTION/MATRIX_APRIORI L COVA\n')

        # Write SOLUTION/MATRIX_APRIORI block (case of normal matrix)
        if ((type(snx.Nc).__name__ != 'NoneType') and (not('matrices' in dont_write)) and (not('apriori' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/MATRIX_APRIORI L INFO\n')
            f.write('*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______\n')
            for i in range(len(snx.prior)):
                j = 0
                while (j <= i):
                    if (snx.Nc[i,j] != 0):
                        f.write(' {0:5} {1:5} {2:21.14e}'.format(i+1, j+1, snx.Nc[i,j]))
                        if ((j+2 <= i) and (snx.Nc[i,j+2] != 0)):
                            f.write(' {0:21.14e}'.format(snx.Nc[i,j+1]))
                            f.write(' {0:21.14e}\n'.format(snx.Nc[i,j+2]))
                            j = j+3
                        elif ((j+1 <= i) and (snx.Nc[i,j+1] != 0)):
                            f.write(' {0:21.14e}\n'.format(snx.Nc[i,j+1]))
                            j = j+2
                        else:
                            f.write('\n')
                            j = j+1
                    else:
                        j = j+1
            f.write('-SOLUTION/MATRIX_APRIORI L INFO\n')

        # Write SOLUTION/MATRIX_ESTIMATE block (case of covariance matrix)
        if ((type(snx.Q).__name__ != 'NoneType') and (not('matrices' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/MATRIX_ESTIMATE L COVA\n')
            f.write('*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______\n')
            for i in range(snx.npar):
                #j = 0
                #while (j <= i):
                    #f.write(' {0:5} {1:5} {2:21.14e}'.format(i+1, j+1, snx.Q[i,j]))
                    #if (j+1 <= i):
                        #j = j+1
                        #f.write(' {0:21.14e}'.format(snx.Q[i,j]))
                        #if (j+1 <= i):
                            #j = j+1
                            #f.write(' {0:21.14e}\n'.format(snx.Q[i,j]))
                        #else:
                            #f.write('\n')
                    #else:
                        #f.write('\n')
                    #j = j+1
                j = 0
                while (j <= i):
                    if (snx.Q[i,j] != 0):
                        f.write(' {0:5} {1:5} {2:21.14e}'.format(i+1, j+1, snx.Q[i,j]))
                        if ((j+2 <= i) and (snx.Q[i,j+2] != 0)):
                            f.write(' {0:21.14e}'.format(snx.Q[i,j+1]))
                            f.write(' {0:21.14e}\n'.format(snx.Q[i,j+2]))
                            j = j+3
                        elif ((j+1 <= i) and (snx.Q[i,j+1] != 0)):
                            f.write(' {0:21.14e}\n'.format(snx.Q[i,j+1]))
                            j = j+2
                        else:
                            f.write('\n')
                            j = j+1
                    else:
                        j = j+1
            f.write('-SOLUTION/MATRIX_ESTIMATE L COVA\n')

        # Write SOLUTION/MATRIX_ESTIMATE block (case of total normal matrix)
        if ((snx.Ntot != None) and (not('matrices' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/MATRIX_ESTIMATE L INFO\n')
            f.write('*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______\n')
            for i in range(snx.npar):
                #j = 0
                #while (j <= i):
                    #f.write(' {0:5} {1:5} {2:21.14e}'.format(i+1, j+1, snx.Ntot[i,j]))
                    #if (j+1 <= i):
                        #j = j+1
                        #f.write(' {0:21.14e}'.format(snx.Ntot[i,j]))
                        #if (j+1 <= i):
                            #j = j+1
                            #f.write(' {0:21.14e}\n'.format(snx.Ntot[i,j]))
                        #else:
                            #f.write('\n')
                    #else:
                        #f.write('\n')
                    #j = j+1
                j = 0
                while (j <= i):
                    if (snx.Ntot[i,j] != 0):
                        f.write(' {0:5} {1:5} {2:21.14e}'.format(i+1, j+1, snx.Ntot[i,j]))
                        if ((j+2 <= i) and (snx.Ntot[i,j+2] != 0)):
                            f.write(' {0:21.14e}'.format(snx.Ntot[i,j+1]))
                            f.write(' {0:21.14e}\n'.format(snx.Ntot[i,j+2]))
                            j = j+3
                        elif ((j+1 <= i) and (snx.Ntot[i,j+1] != 0)):
                            f.write(' {0:21.14e}\n'.format(snx.Ntot[i,j+1]))
                            j = j+2
                        else:
                            f.write('\n')
                            j = j+1
                    else:
                        j = j+1
            f.write('-SOLUTION/MATRIX_ESTIMATE L INFO\n')

        # Write SOLUTION/NORMAL_EQUATION_VECTOR block
        if ((snx.k != None) and (not('matrices' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/NORMAL_EQUATION_VECTOR\n')
            f.write('*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___\n')
            for i in range(snx.npar):
                p =  snx.param[i]
                f.write(' {0:5} {1} {2} {3} {4} '.format(i+1, p.type, p.code, p.pt, p.soln))
                f.write('{0} {1} {2} {3:21.14e}\n'.format(p.tref, p.unit, p.const, snx.k[i]))
            f.write('-SOLUTION/NORMAL_EQUATION_VECTOR\n')

        # Write SOLUTION/NORMAL_EQUATION_MATRIX block
        if ((snx.N != None) and (not('matrices' in dont_write))):
            f.write('*-------------------------------------------------------------------------------\n')
            f.write('+SOLUTION/NORMAL_EQUATION_MATRIX\n')
            f.write('*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______\n')
            for i in range(snx.npar):
                j = 0
                while (j <= i):
                    f.write(' {0:5} {1:5} {2:21.14e}'.format(i+1, j+1, snx.N[i,j]))
                    if (j+1 <= i):
                        j = j+1
                        f.write(' {0:21.14e}'.format(snx.N[i,j]))
                        if (j+1 <= i):
                            j = j+1
                            f.write(' {0:21.14e}\n'.format(snx.N[i,j]))
                        else:
                            f.write('\n')
                    else:
                        f.write('\n')
                    j = j+1
            f.write('-SOLUTION/NORMAL_EQUATION_MATRIX\n')

        # Write last line and close output SINEX file
        f.write('%ENDSNX\n')
        f.close()

#-------------------------------------------------------------------------------
# Routine : check_staid
# Purpose : Check PT codes and DOMES numbers
# Author  : P. Rebischung
# Created : 12-Aug-2011
#
# Changes : PR, 27-Jan-2017 : Add possibility to check station coordinates in case
#                             of multiple possible DOMES numbers
#           PR, 27-Jan-2017 : Add possibility to update SINEX station description
#                             from DOMES number catalogue
#
# Input   : - codomes    : DOMES number catalogue
#           - check_pt   : True if PT codes should be checked. Default is True.
#           - check_crd  : True if station coordinates should be checked in case of
#                          multiple DOMES numbers. Default is True.
#           - check_desc : True if station descriptions should be updated from DOMES
#                          number catalogue. Default is True.
#           - log        : Log file. Default is 'None'.
#           - append     : If true, text will be appended to log file. Default if False.
# Output  : - modif      : True if any modification has to be made
#-------------------------------------------------------------------------------
    def check_staid(self, codomes, check_pt=True, check_crd=True, check_desc=True, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        # If necessary, open log file and write header
        if (log):
            if (append):
                fl = open(log, 'a')
            else:
                fl = open(log, 'w')
            fl.write('================================================================================\n')
            fl.write('sinex::check_staid : Check PT codes and DOMES numbers\n')
            fl.write('================================================================================\n')
            fl.write('\n')

        # Initialization
        modif = False

        # Loop over stations
        for i in range(len(snx.sta)):
            code  = snx.sta[i].code
            pt    = snx.sta[i].pt
            domes = snx.sta[i].domes

            # Check PT code
            # -------------
            if (check_pt):
                mod = False

                # PT should be ' A' for all stations...
                if ((code != 'IISC') and (code != 'KELY') and (pt != ' A')):
                    modif = True
                    mod = True
                    pt2 = ' A'

                # ...except IISC and KELY.
                elif ((code in ['IISC', 'KELY']) and (pt != ' B')):
                    modif = True
                    mod = True
                    pt2 = ' B'

                # If a correction is needed
                if (mod):

                    # Correct PT in snx.sta
                    snx.sta[i].pt = pt2

                    # Correct PT in snx.param
                    for j in range(snx.npar):
                        if ((snx.param[j].code == code) and (snx.param[j].pt == pt)):
                            snx.param[j].pt = pt2

                    # Correct PT in snx.prior
                    if (snx.prior):
                        for j in range(len(snx.prior)):
                            if ((snx.prior[j].code == code) and (snx.prior[j].pt == pt)):
                                snx.prior[j].pt = pt2

                    # Print message in log file
                    if (log):
                        fl.write('PT corrected    : {0} {1} {3} => {0} {2} {3}\n'.format(code, pt, pt2, domes))

            # Check DOMES number
            # ------------------
            mod = False

            # 1st case: Do not check station coordinates
            if (not(check_crd)):
                pass

                ## If code+domes is found in DOMES number catalogue, everything is fine.
                ## Else...
                #if (not(code+domes in [c.code+c.domes for c in codomes])):

                    ## If code is nevertheless found in DOMES number catalogue, a correction is needed.
                    #if (code in [c.code for c in codomes]):
                        #j = [c.code for c in codomes].index(code)
                        #modif = True
                        #mod = True
                        #domes2 = codomes[j].domes

                    ## If code is not found in DOMES number catalogue, set default DOMES number.
                    #elif (domes != default_domes):
                        #modif = True
                        #mod = True
                        #domes2 = default_domes

            # 2nd case: Check station coordinates
            else:

                # Look for every occurence of code in DOMES number catalogue
                ind = numpy.nonzero(numpy.array([c.code for c in codomes]) == code)[0]

                # If no occurence is found, set default DOMES number if needed.
                if (len(ind) == 0):
                    if (domes != default_domes):
                        modif = True
                        mod = True
                        domes2 = default_domes

                # If one occurence is found,
                elif (len(ind) == 1):
                    ind = ind[0]

                    # Change DOMES number if needed
                    if (domes != codomes[ind].domes):
                        modif = True
                        mod = True
                        domes2 = codomes[ind].domes

                    # Change station description if needed
                    if (check_desc):
                        snx.sta[i].description = codomes[ind].description

                # If several occurences are found,
                else:

                    # Get station coordinates
                    (lon, lat) = snx.get_lonlat([code])
                    lam = pi/180*lon[0]
                    phi = pi/180*lat[0]

                    # Compute distances to every point of the DOMES number catalogue
                    d = numpy.zeros(len(ind))
                    for j in range(len(ind)):
                        lamj = pi/180*codomes[ind[j]].lon
                        phij = pi/180*codomes[ind[j]].lat
                        d[j] = acos(sin(phi)*sin(phij) + cos(phi)*cos(phij)*cos(lam-lamj))

                    # Index of the closest point
                    j = numpy.nonzero(d == numpy.min(d))[0][0]

                    # Change DOMES number if needed
                    if (domes != codomes[ind[j]].domes):
                        modif = True
                        mod = True
                        domes2 = codomes[ind[j]].domes

                    # Change station description if needed
                    if (check_desc):
                        snx.sta[i].description = codomes[ind[j]].description

            # Correct DOMES number in snx.sta
            if (mod):
                snx.sta[i].domes = domes2

                # Print message in log file
                if (log):
                    fl.write('DOMES corrected : {0} {1} {2} => {0} {1} {3}\n'.format(code, pt, domes, domes2))

        # Close log file if necessary
        if (log):
            fl.write('\n')
            fl.close()

        return modif

#-------------------------------------------------------------------------------
# Routine : check_solns
# Purpose : Check solns in an "instantaneous" solution
# Author  : P. Rebischung
# Created : 12-Aug-2011
#
# Changes :
#
# Input   : - solns  : Discontinuity list
#           - log    : Log file. Default is 'None'.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def check_solns(self, solns, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        # If necessary, open log file and write header
        if (log):
            if (append):
                fl = open(log, 'a')
            else:
                fl = open(log, 'w')
            fl.write('================================================================================\n')
            fl.write('sinex::check_solns : Check solution numbers\n')
            fl.write('================================================================================\n')
            fl.write('\n')

        # 1 - Check solns of estimated parameters
        # ---------------------------------------
        for i in range(snx.npar):
            if (snx.param[i].type == 'STAX  '):
                p = snx.param[i]

                # If current station is found in solns table,
                if (p.code+p.pt in [s.code+s.pt for s in solns]):
                    ista = [s.code+s.pt for s in solns].index(p.code+p.pt)

                    # Look for appropriate soln
                    t = snx.param[i].tref
                    isoln = 0
                    while ((solns[ista].P[isoln].end != '00:000:00000') and (earlier(solns[ista].P[isoln].end, t))):
                        isoln = isoln+1
                    soln2 = solns[ista].P[isoln].soln

                # Else, default soln is '   1'
                else:
                    soln2 = '   1'
                    if (log):
                        fl.write('### {0.code} {0.pt} not found in soln catalogue !\n'.format(p))

                # If soln has to be modified,
                if (snx.param[i].soln != soln2):

                    # Print message in log file
                    if (log):
                        fl.write('{0.code} {0.pt} {0.soln} => {0.code} {0.pt} {1}\n'.format(p, soln2))

                    # Modify soln in snx.param
                    snx.param[i+0].soln = soln2
                    snx.param[i+1].soln = soln2
                    snx.param[i+2].soln = soln2

                # Modify soln in snx.sta
                ista = [s.code+s.pt for s in snx.sta].index(p.code+p.pt)
                snx.sta[ista].soln[0].soln = soln2

        # 2 - Check solns of a priori parameters
        # --------------------------------------
        if (snx.prior):
            for i in range(len(snx.prior)):
                if (snx.prior[i].type == 'STAX  '):
                    p = snx.prior[i]

                    # If current station is found in solns table,
                    if (p.code+p.pt in [s.code+s.pt for s in solns]):
                        ista = [s.code+s.pt for s in solns].index(p.code+p.pt)

                        # Look for appropriate soln
                        t = snx.prior[i].tref
                        isoln = 0
                        while ((solns[ista].P[isoln].end != '00:000:00000') and (earlier(solns[ista].P[isoln].end, t))):
                            isoln = isoln+1
                        soln2 = solns[ista].P[isoln].soln

                    # Else, default soln is '   1'
                    else:
                        soln2 = '   1'

                    # If necessary, modify soln in snx.prior
                    if (snx.prior[i].soln != soln2):
                        snx.prior[i+0].soln = soln2
                        snx.prior[i+1].soln = soln2
                        snx.prior[i+2].soln = soln2

        # Close log file if necessary
        if (log):
            fl.write('\n')
            fl.close()

#-------------------------------------------------------------------------------
# Routine : check_erp
# Purpose : Set ERP reference epochs to exactly noon
# Author  : P. Rebischung
# Created : 24-May-2012
#
# Changes :
#
# Input   :
# Output  :
#-------------------------------------------------------------------------------
    def check_erp(self):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        # Check epochs of estimated ERPs
        for p in snx.param:
            if (p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']):
                if (p.tref[-5:] != '43200'):
                    p.tref = p.tref[:-5] + '43200'

        # Check epochs of a priori ERPs
        if (snx.prior):
            for p in snx.prior:
                if (p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']):
                    if (p.tref[-5:] != '43200'):
                        p.tref = p.tref[:-5] + '43200'

#-------------------------------------------------------------------------------
# Routine : get_xyz
# Purpose : Get cartesian coordinates of specified stations
# Author  : P. Rebischung
# Created : 05-Jul-2011
#
# Changes :
#
# Input   : - code : 4-char codes
#           - pt   : PT codes. Default is None.
#           - soln : Solns. Default is None.
# Output  : - X    : [X, Y, Z] cartesian coordinates in m
#-------------------------------------------------------------------------------
    def get_xyz(self, code, pt=None, soln=None):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        # Initialization
        X = numpy.zeros((len(code), 3))

        # Loop over requested points
        for i in range(len(code)):

            # Look for a STAX parameter for current point
            f = 0
            j = 0
            if ((not(pt)) and (not(soln))):
                if ('STAX  '+code[i] in [p.type+p.code for p in snx.param]):
                    f = 1
                    j = [p.type+p.code for p in snx.param].index('STAX  '+code[i])
            elif (not(soln)):
                if ('STAX  '+code[i]+pt[i] in [p.type+p.code+p.pt for p in snx.param]):
                    f = 1
                    j = [p.type+p.code+p.pt for p in snx.param].index('STAX  '+code[i]+pt[i])
            else:
                if ('STAX  '+code[i]+pt[i]+soln[i] in [p.type+p.code+p.pt+p.soln for p in snx.param]):
                    f = 1
                    j = [p.type+p.code+p.pt+p.soln for p in snx.param].index('STAX  '+code[i]+pt[i]+soln[i])

            # If a STAX parameter was found for current point, add its coordinates into X
            if (f):
                X[i] = [snx.x[j+0], snx.x[j+1], snx.x[j+2]]

        return X

#-------------------------------------------------------------------------------
# Routine : get_plh
# Purpose : Get geographical coordinates of specified stations
# Author  : P. Rebischung
# Created : 05-Jul-2011
#
# Changes :
#
# Input   : - code : 4-char codes
#           - pt   : PT codes. Default is None.
#           - soln : Solns. Default is None.
# Output  : - phi  : Latitudes (rad)
#           - lam  : Longitudes (rad)
#           - h    : Ellipsoidal heights (m)
#-------------------------------------------------------------------------------
    def get_plh(self, code, pt=None, soln=None):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        X = snx.get_xyz(code, pt, soln)

        return cart2geo(X)

#-------------------------------------------------------------------------------
# Routine : get_lonlat
# Purpose : Get longitudes and latitudes of specified stations
# Author  : P. Rebischung
# Created : 05-Jul-2011
#
# Changes :
#
# Input   : - code : 4-char codes
#           - pt   : PT codes. Default is None.
#           - soln : Solns. Default is None.
# Output  : - lon  : Longitudes (deg)
#           - lat  : Latitudes (deg)
#-------------------------------------------------------------------------------
    def get_lonlat(self, code, pt=None, soln=None):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        (phi, lam, h) = snx.get_plh(code, pt, soln)

        return (180/pi*lam, 180/pi*phi)

#-------------------------------------------------------------------------------
# Routine : get_sigenh
# Purpose : Get E,N,H formal errors of specified stations
# Author  : P. Rebischung
# Created : 20-May-2012
#
# Changes :
#
# Input   : - code : 4-char codes
#           - pt   : PT codes
#           - soln : Solns
# Output  :
#-------------------------------------------------------------------------------
    def get_sigenh(self, code, pt, soln):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Station positions
        X = snx.get_xyz(code, pt, soln)

        # Geocentric to topocentric rotation matrices
        R = xyz2enh(X)

        # Get E,N,H formal errors
        s = numpy.zeros((len(code), 3))
        for i in range(len(code)):
            ind = [p.type+p.code+p.pt+p.soln for p in snx.param].index('STAX  '+code[i]+pt[i]+soln[i])
            Q = snx.Q[ind:ind+3][:,ind:ind+3]
            Q = dot(dot(R[i], Q), R[i].T)
            s[i,:] = numpy.sqrt(numpy.diag(Q))

        return s


#-------------------------------------------------------------------------------
# Routine : get_cov_sta
# Purpose : Get the site covariance matrix for a list of sites
# Author  : J.-M. Nocquet
# Created : 26-September-2018
#
# Changes :
#
# Input   : - code : 4-char codes
#           - pt   : PT codes
#           - soln : Solns
# Output  :
#-------------------------------------------------------------------------------
    def get_cov_sta(self, code, pt=None, soln=None):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        lQ = []

        # Get individual site covariance from the overall solution covariance matrix

        for i in range(len(code)):

            if ( (pt is None ) and (soln is None) ):
                if ('STAX  '+code[i] in [p.type+p.code for p in snx.param]):
                    ind = [p.type+p.code for p in snx.param].index('STAX  '+code[i])
            
            lQ.append( snx.Q[ind:ind+3][:,ind:ind+3] )

        return lQ

#-------------------------------------------------------------------------------
# Routine : get_helmert_matrix
# Purpose : Get design matrix of Helmert parameters
# Author  : P. Rebischung
# Created : 02-Sep-2014
#
# Changes :
#
# Input   : - params : String indicating which Helmert parameters should be considered.
#                      It can be composed by any combination of letters 'T' (translations),
#                     'S' (scale) and 'R' (rotations). Default is 'RST' (all 7 parameters).
#           - pole   : If True, the RY/XPO and RX/YPO partial derivatives are set to 1.
#                      Default is False.
#           - gc     : If True, the TX/XGC, TY/YGC and TZ/ZGC partial derivatives are set to 1.
# Output  : - A      : Design matrix of Helmert parameters
#-------------------------------------------------------------------------------
    def get_helmert_matrix(self, params='RST', pole=False, gc=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initialization
        A = numpy.zeros((snx.npar, 7))

        # Loop over STAX parameters in order to fill the full design matrix A
        for i in range(snx.npar):
            p = snx.param[i]
            if (p.type == 'STAX  '):
                A[i+0, 0]  =  1e-3
                A[i+1, 1]  =  1e-3
                A[i+2, 2]  =  1e-3
                A[i+0, 3]  =  1e-9 * snx.x[i+0]
                A[i+1, 3]  =  1e-9 * snx.x[i+1]
                A[i+2, 3]  =  1e-9 * snx.x[i+2]
                A[i+1, 4]  = -mas2rad * snx.x[i+2]
                A[i+2, 4]  =  mas2rad * snx.x[i+1]
                A[i+0, 5]  =  mas2rad * snx.x[i+2]
                A[i+2, 5]  = -mas2rad * snx.x[i+0]
                A[i+0, 6]  = -mas2rad * snx.x[i+1]
                A[i+1, 6]  =  mas2rad * snx.x[i+0]
            elif ((p.type == 'XPO   ') and (pole)):
                A[i, 5] = 1.
            elif ((p.type == 'YPO   ') and (pole)):
                A[i, 4] = 1.
            elif ((p.type == 'XGC   ') and (gc)):
                A[i, 0] = 1e-3
            elif ((p.type == 'YGC   ') and (gc)):
                A[i, 1] = 1e-3
            elif ((p.type == 'ZGC   ') and (gc)):
                A[i, 2] = 1e-3

        # Keep relevant columns of A
        ind = []
        if ('T' in params):
            ind.extend(list(range(0, 3)))
        if ('S' in params):
            ind.append(3)
        if ('R' in params):
            ind.extend(list(range(4, 7)))
        A = A[:,ind]

        return A

#-------------------------------------------------------------------------------
# Routine : map
# Purpose : Draw station map
# Author  : P. Rebischung
# Created : 08-Feb-2012
#
# Changes :
#
# Input   : - write_codes : True for station codes to be written on the map. Default is True.
#           - title       : Map title. Default is None.
#           - output      : Output file. Default is None.
# Output  :
#-------------------------------------------------------------------------------
    def map(self, write_codes=True, title=None, output=None):


        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        code = [s.code for s in snx.sta]
        (lon, lat) = snx.get_lonlat(code)
        station_map(lon, lat, code, write_codes, title, output)

#-------------------------------------------------------------------------------
# Routine : propagate
# Purpose : Propagate station positions to specified date
# Author  : P. Rebischung
# Created : 23-May-2011
#
# Changes :
#
# Input   : - epo      : Epoch (SINEX format)
#           - keep_vel : True or False depending on whether station velocities
#                        should be kept or not. Default is False.
#-------------------------------------------------------------------------------
    def propagate(self, epo, keep_vel=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # date object at propagation epoch
        t = date.from_snxepoch(epo)

        # Loop over parameters to
        # - build the design matrix A
        # - build a table ind containing the indices of non-velocity parameters
        A = sparse.identity(snx.npar, format='lil')
        ind = []
        for i in range(snx.npar):

            # If parameter i is a velocity,
            if (snx.param[i].type[0:3] == 'VEL'):

                # Current reference epoch of corresponding 'STA' parameter => dt in years
                ti = date.from_snxepoch(snx.param[i-3].tref)
                dti = (t.mjd - ti.mjd) / 365.25

                # Fill design matrix
                A[i-3, i] = dti

            # Else, add parameter i to the list of non-velocity parameters.
            elif (snx.param[i].type[0:3] == 'STA'):
                ind.append(i)

        # Loop over STA and VEL parameters to change their reference dates
        for p in snx.param:
            if (p.type[0:3] in ['STA', 'VEL']):
                p.tref = epo

        # If velocities should not be kept in sinex object, make some cleaning.
        A = A.tocsr()
        if (not(keep_vel)):
            A = A[ind,:]
            snx.npar = len(ind)
            snx.param = [snx.param[i] for i in ind]

        # Propagate coordinates
        snx.x = A.dot(snx.x)

        # Propagate covariance matrix if available
        if (snx.Q != None):
            snx.Q = A.dot((A.dot(snx.Q)).T)

        # Update standard deviations if covariance matrix is available
        # Else, set them to 0.
        if (snx.Q != None):
            snx.sig = numpy.sqrt(numpy.diag(snx.Q))
        else:
            snx.sig = numpy.zeros(snx.npar)

#-------------------------------------------------------------------------------
# Routine : remove_sta
# Purpose : Remove specified stations from a solution
# Author  : P. Rebischung
# Created : 16-Sep-2011
#
# Changes :
#
# Input   : - code   : 4-char codes
#           - pt     : PT codes. Default is None.
#           - soln   : Solns. Default is None.
#           - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def remove_sta(self, code, pt=None, soln=None, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # If necessary, open log file and write header
        if (log):
            if (append):
                fl = open(log, 'a')
            else:
                fl = open(log, 'w')
            fl.write('================================================================================\n')
            fl.write('sinex::remove_sta : Remove specified stations from a solution\n')
            fl.write('================================================================================\n')
            fl.write('\n')

        # 1 - Deal with estimated parameters and snx.sta
        #-----------------------------------------------

        # Initialization
        indk = []

        # Loop over parameters
        for i in range(snx.npar):
            p = snx.param[i]

            # Must parameter i be kept?
            keep = False
            if ((pt) and (soln)):
                if (not(p.code+p.pt+p.soln in [code[j]+pt[j]+soln[j] for j in range(len(code))])):
                    keep = True
                    indk.append(i)
            elif (pt):
                if (not(p.code+p.pt in [code[j]+pt[j] for j in range(len(code))])):
                    keep = True
                    indk.append(i)
            else:
                if (not(p.code in code)):
                    keep = True
                    indk.append(i)

            # If parameter i has to be removed, print message in log file.
            if ((not(keep)) and (log)):
                fl.write('{0.type} {0.code} {0.pt} {0.soln} removed\n'.format(p))

        # If any parameter has to be deleted,
        if (len(indk) < snx.npar):

            # Update snx.npar, snx.param, snx.x, snx.sig and snx.Q
            snx.npar  = len(indk)
            snx.param = [snx.param[i] for i in indk]
            snx.x     = snx.x[indk]
            snx.sig   = snx.sig[indk]
            if (snx.Q != None):
                snx.Q = snx.Q[numpy.ix_(indk,indk)]

            # List of points for which there remain estimated positions
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'STAX  ')[0]
            codek = [snx.param[i].code for i in ind]
            ptk = [snx.param[i].pt for i in ind]
            solnk = [snx.param[i].soln for i in ind]

            # Update snx.sta[*].soln
            for s in snx.sta:
                i = 0
                while (i < len(s.soln)):
                    if (not(s.code+s.pt+s.soln[i].soln in [codek[k]+ptk[k]+solnk[k] for k in range(len(codek))])):
                        (s.soln).pop(i)
                    else:
                        i = i+1

            # Update snx.sta
            i = 0
            while (i < len(snx.sta)):
                if (len(snx.sta[i].soln) == 0):
                    (snx.sta).pop(i)
                else:
                    i = i+1

            ## PR, 02-Feb-2015 : Simplified handling of snx.sta
            #i = 0
            #while (i < len(snx.sta)):
                #if (not(snx.sta[i].code+snx.sta[i].pt in [codek[k]+ptk[k] for k in range(len(codek))])):
                    #(snx.sta).pop(i)
                #else:
                    #i = i+1

        # 2 - Deal with a priori parameters
        #----------------------------------

        if (snx.prior):

            # Initialization
            indk = []

            # Loop over a priori parameters
            for i in range(len(snx.prior)):
                p = snx.prior[i]

                # Must a priori parameter i be kept?
                keep = False
                if ((pt) and (soln)):
                    if (not(p.code+p.pt+p.soln in [code[j]+pt[j]+soln[j] for j in range(len(code))])):
                        keep = True
                        indk.append(i)
                elif (pt):
                    if (not(p.code+p.pt in [code[j]+pt[j] for j in range(len(code))])):
                        keep = True
                        indk.append(i)
                else:
                    if (not(p.code in code)):
                        keep = True
                        indk.append(i)

            # If any parameter has to be deleted,
            if (len(indk) < len(snx.prior)):

                # Update snx.prior, snx.x0, snx.sig0, snx.Qc and snx.Nc
                snx.prior = [snx.prior[i] for i in indk]
                snx.x0    = snx.x0[indk]
                snx.sig0  = snx.sig0[indk]
                if (snx.Qc != None):
                    snx.Qc = snx.Qc[numpy.ix_(indk,indk)]
                if (snx.Nc != None):
                    snx.Nc = snx.Nc[numpy.ix_(indk,indk)]

        # Close log file if necessary
        if (log):
            fl.write('\n')
            fl.close()

#-------------------------------------------------------------------------------
# Routine : keep_sta
# Purpose : Keep only specified stations in a solution
# Author  : P. Rebischung
# Created : 16-Sep-2011
#
# Changes : PR, 05-Jul-2016 : Deal with PSD parameters
#
# Input   : - code   : 4-char codes
#           - pt     : PT codes. Default is None.
#           - soln   : Solns. Default is None.
#           - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def keep_sta(self, code, pt=None, soln=None, log=None, append=False):


        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        # If necessary, open log file and write header
        if (log):
            if (append):
                fl = open(log, 'a')
            else:
                fl = open(log, 'w')
            fl.write('================================================================================\n')
            fl.write('sinex::keep_sta : Keep only specified stations in a solution\n')
            fl.write('================================================================================\n')
            fl.write('\n')

        # 1 - Deal with estimated parameters and snx.sta
        #-----------------------------------------------

        # Initialization
        indk = []

        # Loop over parameters
        for i in range(snx.npar):
            p = snx.param[i]

            # Must parameter i be kept?
            keep = False
            if ((not(p.type[0:3] in ['STA', 'VEL'])) and (not(p.type[0:4] in ['AEXP', 'TEXP', 'ALOG', 'TLOG']))):
                keep = True
                indk.append(i)
            elif ((pt) and (soln)):
                if (p.code+p.pt+p.soln in [code[j]+pt[j]+soln[j] for j in range(len(code))]):
                    keep = True
                    indk.append(i)
            elif (pt):
                if (p.code+p.pt in [code[j]+pt[j] for j in range(len(code))]):
                    keep = True
                    indk.append(i)
            else:
                if (p.code in code):
                    keep = True
                    indk.append(i)

            # If parameter i has to be removed, print message in log file.
            if ((not(keep)) and (log)):
                fl.write('{0.type} {0.code} {0.pt} {0.soln} removed\n'.format(p))

        # If any parameter has to be deleted,
        if (len(indk) < snx.npar):

            # Update snx.npar, snx.param, snx.x, snx.sig and snx.Q
            snx.npar  = len(indk)
            snx.param = [snx.param[i] for i in indk]
            snx.x     = snx.x[indk]
            snx.sig   = snx.sig[indk]
            if (snx.Q != None):
                snx.Q = snx.Q[numpy.ix_(indk,indk)]

            # List of points for which there remain estimated positions
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'STAX  ')[0]
            codek = [snx.param[i].code for i in ind]
            ptk = [snx.param[i].pt for i in ind]
            solnk = [snx.param[i].soln for i in ind]

            # Update snx.sta[*].soln
            for s in snx.sta:
                i = 0
                while (i < len(s.soln)):
                    if (not(s.code+s.pt+s.soln[i].soln in [codek[k]+ptk[k]+solnk[k] for k in range(len(codek))])):
                        (s.soln).pop(i)
                    else:
                        i = i+1

            # Update snx.sta
            i = 0
            while (i < len(snx.sta)):
                if (len(snx.sta[i].soln) == 0):
                    (snx.sta).pop(i)
                else:
                    i = i+1

            ## PR, 02-Feb-2015 : Simplified handling of snx.sta
            #i = 0
            #while (i < len(snx.sta)):
                #if (not(snx.sta[i].code+snx.sta[i].pt in [code[k]+pt[k] for k in range(len(code))])):
                    #(snx.sta).pop(i)
                #else:
                    #i = i+1

        # 2 - Deal with a priori parameters
        #----------------------------------

        if (snx.prior):

            # Initialization
            indk = []

            # Loop over a priori parameters
            for i in range(len(snx.prior)):
                p = snx.prior[i]

                # Must parameter i be kept?
                keep = False
                if ((not(p.type[0:3] in ['STA', 'VEL'])) and (not(p.type[0:4] in ['AEXP', 'TEXP', 'ALOG', 'TLOG']))):
                    keep = True
                    indk.append(i)
                elif ((pt) and (soln)):
                    if (p.code+p.pt+p.soln in [code[j]+pt[j]+soln[j] for j in range(len(code))]):
                        keep = True
                        indk.append(i)
                elif (pt):
                    if (p.code+p.pt in [code[j]+pt[j] for j in range(len(code))]):
                        keep = True
                        indk.append(i)
                else:
                    if (p.code in code):
                        keep = True
                        indk.append(i)

            # If any parameter has to be deleted,
            if (len(indk) < len(snx.prior)):

                # Update snx.prior, snx.x0, snx.sig0, snx.Qc and snx.Nc
                snx.prior = [snx.prior[i] for i in indk]
                snx.x0    = snx.x0[indk]
                snx.sig0  = snx.sig0[indk]
                if (snx.Qc != None):
                    snx.Qc = snx.Qc[numpy.ix_(indk,indk)]
                if (snx.Nc != None):
                    snx.Nc = snx.Nc[numpy.ix_(indk,indk)]

        # Close log file if necessary
        if (log):
            fl.write('\n')
            fl.close()

#-------------------------------------------------------------------------------
# Routine : remove_sta_wo_domes
# Purpose : Remove stations without DOMES numbers from a solution
# Author  : P. Rebischung
# Created : 19-May-2012
#
# Changes :
#
# Input   :
# Output  :
#-------------------------------------------------------------------------------
    def remove_sta_wo_domes(self):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initializations
        code = []
        pt   = []

        # Get list of stations without DOMES numbers
        for s in snx.sta:
            if (s.domes == default_domes):
                code.append(s.code)
                pt.append(s.pt)

        # Remove those stations
        snx.remove_sta(code, pt)

#-------------------------------------------------------------------------------
# Routine : keep_relevant_solns
# Purpose : Extract solns that a relevant to specified epoch from a solution
# Author  : P. Rebischung
# Created : 27-Jun-2012
#
# Changes :
#
# Input   : - t     : Epoch (SINEX format)
#           - solns : Soln table
# Output  :
#-------------------------------------------------------------------------------
    def keep_relevant_solns(self, t, solns):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initializations
        code = []
        pt   = []
        soln = []

        # Loop over STAX parameters
        for p in snx.param:
            if (p.type == 'STAX  '):

                # If current soln is in reference soln table,
                if (p.code+p.pt in [r.code+r.pt for r in solns]):
                    i = [r.code+r.pt for r in solns].index(p.code+p.pt)
                    if (p.soln in [r.soln for r in solns[i].P]):
                        j = [r.soln for r in solns[i].P].index(p.soln)
                        start = solns[i].P[j].start
                        end   = solns[i].P[j].end

                        # And if it is relevant to requested epoch,
                        if (((start == '00:000:00000') or (earlier(start, t))) and ((end == '00:000:00000') or (earlier(t, end)))):

                            # Add it to the list of solns to keep
                            code.append(p.code)
                            pt.append(p.pt)
                            soln.append(p.soln)

        # Extract relevant solns
        snx.keep_sta(code, pt, soln)

#-------------------------------------------------------------------------------
# Routine : remove_params
# Purpose : Remove certain types of parameters from a solution
# Author  : P. Rebischung
# Created : 11-Aug-2011
#
# Changes :
#
# Input   : types : List of keywords to indicate on which parameters constraints should be removed.
#                   This can include the following keywords:
#                   'XPO', 'XPOR', 'YPO', 'YPOR', 'UT', 'LOD'; 'ERP' for all kinds of ERPs;
#                   'SATA_X', 'SATA_Y', 'SATA_Z'; 'SATA' for all satellite parameters;
#                   'STA'; 'VEL'; 'GC' and 'TRANS' for all transformation parameters.
# Output  :
#-------------------------------------------------------------------------------
    def remove_params(self, types):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # 1 - Deal with estimated parameters
        #-----------------------------------

        # Get indices of parameters to keep
        indk = []
        for i in range(snx.npar):
            p = snx.param[i]
            if not(((p.type == 'XPO   ') and ('XPO'  in types))  or
                   ((p.type == 'XPOR  ') and ('XPOR'  in types)) or
                   ((p.type == 'YPO   ') and ('YPO'  in types))  or
                   ((p.type == 'YPOR  ') and ('YPOR'  in types)) or
                   ((p.type == 'UT    ') and ('UT'  in types))   or
                   ((p.type == 'LOD   ') and ('LOD'  in types))  or
                   ((p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']) and ('ERP' in types)) or
                   ((p.type == 'SATA_X') and ('SATA_X' in types)) or
                   ((p.type == 'SATA_Y') and ('SATA_Y' in types)) or
                   ((p.type == 'SATA_Z') and ('SATA_Z' in types)) or
                   ((p.type[0:4] == 'SATA') and ('SATA' in types)) or
                   ((p.type[0:3] == 'STA') and ('STA' in types)) or
                   ((p.type[0:3] == 'VEL') and ('VEL' in types)) or
                   ((p.type[1:3] == 'GC') and ('GC' in types)) or
                   ((p.type in ['TX    ', 'TY    ', 'TZ    ', 'SC    ', 'RX    ', 'RY    ', 'RZ    ']) and ('TRANS' in types))):
                indk.append(i)

        # Update sinex object
        if (len(indk) < snx.npar):
            snx.npar  = len(indk)
            snx.param = [snx.param[i] for i in indk]
            snx.x     = snx.x[indk]
            snx.sig   = snx.sig[indk]
            if (snx.Q != None):
                snx.Q = snx.Q[numpy.ix_(indk,indk)]
            if ('STA' in types):
                snx.sta = None

        # 2 - Deal with a priori parameters
        #----------------------------------

        if (snx.prior):

            # Get indices of parameters to keep
            indk = []
            for i in range(len(snx.prior)):
                p = snx.prior[i]
                if not(((p.type == 'XPO   ') and ('XPO'  in types))  or
                      ((p.type == 'XPOR  ') and ('XPOR'  in types)) or
                      ((p.type == 'YPO   ') and ('YPO'  in types))  or
                      ((p.type == 'YPOR  ') and ('YPOR'  in types)) or
                      ((p.type == 'UT    ') and ('UT'  in types))   or
                      ((p.type == 'LOD   ') and ('LOD'  in types))  or
                      ((p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']) and ('ERP' in types)) or
                      ((p.type == 'SATA_X') and ('SATA_X' in types)) or
                      ((p.type == 'SATA_Y') and ('SATA_Y' in types)) or
                      ((p.type == 'SATA_Z') and ('SATA_Z' in types)) or
                      ((p.type[0:4] == 'SATA') and ('SATA' in types)) or
                      ((p.type[0:3] == 'STA') and ('STA' in types)) or
                      ((p.type[0:3] == 'VEL') and ('VEL' in types)) or
                      ((p.type[1:3] == 'GC') and ('GC' in types)) or
                      ((p.type in ['TX    ', 'TY    ', 'TZ    ', 'SC    ', 'RX    ', 'RY    ', 'RZ    ']) and ('TRANS' in types))):
                    indk.append(i)

            # Update sinex object
            if (len(indk) < len(snx.prior)):
                snx.prior = [snx.prior[i] for i in indk]
                snx.x0    = snx.x0[indk]
                snx.sig0  = snx.sig0[indk]
                if (snx.Qc != None):
                    snx.Qc = snx.Qc[numpy.ix_(indk,indk)]
                if (snx.Nc != None):
                    snx.Nc = snx.Nc[numpy.ix_(indk,indk)]

#-------------------------------------------------------------------------------
# Routine : remove_undef_params
# Purpose : Remove "undefined" parameters from a solution
# Author  : P. Rebischung
# Created : 11-Aug-2011
#
# Changes :
#
# Input   :
# Output  :
#-------------------------------------------------------------------------------
    def remove_undef_params(self):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Return if sinex object doesn't contain a priori information
        if (not(snx.prior)):
            return

        # Initializations
        iestdel = []
        iaprdel = []

        # Loop over estimated parameters
        for i in range(snx.npar):
            p = snx.param[i]

            # Look for current parameter in snx.prior
            j = -1
            if (p.type+p.code+p.pt+p.soln in [p0.type+p0.code+p0.pt+p0.soln for p0 in snx.prior]):
                j = [p0.type+p0.code+p0.pt+p0.soln for p0 in snx.prior].index(p.type+p.code+p.pt+p.soln)

            # If current parameter has to be deleted,
            if ((snx.sig[i] == 0) or ((snx.param[i].type[0:4] != 'SATA') and (j >= 0) and (snx.sig0[j] > 0) and (snx.sig[i] == snx.sig0[j]))):
                iestdel.append(i)
                if (j >= 0):
                    iaprdel.append(j)

        # Update estimated parameters
        if (len(iestdel) > 0):
            indk = numpy.setdiff1d(list(range(snx.npar)), iestdel)
            snx.npar  = len(indk)
            snx.param = [snx.param[i] for i in indk]
            snx.x     = snx.x[indk]
            snx.sig   = snx.sig[indk]
            if (snx.Q != None):
                snx.Q = snx.Q[numpy.ix_(indk,indk)]

        # Update a priori parameters
        if (len(iaprdel) > 0):
            indk = numpy.setdiff1d(list(range(len(snx.prior))), iaprdel)
            snx.prior = [snx.prior[i] for i in indk]
            snx.x0    = snx.x0[indk]
            snx.sig0  = snx.sig0[indk]
            if (snx.Qc != None):
                snx.Qc = snx.Qc[numpy.ix_(indk,indk)]
            elif (snx.Nc != None):
                snx.Nc = snx.Nc[numpy.ix_(indk,indk)]

#-------------------------------------------------------------------------------
# Routine : unconstrain
# Purpose : Recover unconstrained normal equation from a solution
# Author  : P. Rebischung
# Created : 11-Aug-2011
#
# Changes :
#
# Input   :
# Output  :
#-------------------------------------------------------------------------------
    def unconstrain(self):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # If no a priori information is available, just invert covariance matrix.
        if (not(snx.prior)):
            snx.prior = copy.deepcopy(snx.param)
            snx.x0    = snx.x.copy()
            snx.sig0  = numpy.zeros(snx.npar)
            snx.k     = numpy.zeros(snx.npar)
            snx.N     = syminv(snx.Q)
            snx.Q     = None
            return

        # Get indices of a priori parameters which really correspond to estimated parameters (ind2).
        # Also get indices of corresponding estimated parameters (ind1).
        ind1 = []
        ind2 = []
        for i in range(len(snx.prior)):
            p0 = snx.prior[i]

            # Look for estimated parameter corresponding to snx.prior[i]
            if (p0.type+p0.code+p0.pt+p0.soln in [p.type+p.code+p.pt+p.soln for p in snx.param]):
                j = [p.type+p.code+p.pt+p.soln for p in snx.param].index(p0.type+p0.code+p0.pt+p0.soln)
                ind1.append(j)
                ind2.append(i)

        # If no a priori parameter corresponds to an estimated parameter, just invert covariance matrix.
        if (len(ind2) == 0):
            snx.prior = copy.deepcopy(snx.param)
            snx.x0    = snx.x.copy()
            snx.sig0  = numpy.zeros(snx.npar)
            snx.k     = numpy.zeros(snx.npar)
            snx.N     = syminv(snx.Q)
            snx.Q     = None
            snx.Qc    = None
            snx.Nc    = None
            return

        # Get normal matrix of constraints
        Nc = None
        if (type(snx.Nc).__name__ != 'NoneType'):
            Nc = snx.Nc
            snx.Nc = None
        elif (type(snx.Qc).__name__ != 'NoneType'):
            Nc = sympinv(snx.Qc)
            snx.Qc = None
        else:
            Nc = numpy.zeros((len(snx.prior), len(snx.prior)))

        # Complete diagonal of the normal matrix of constraints with 1/sig0**2
        # for potential a priori parameters which do not appear in the MATRIX_APRIORI block
        for i in range(len(snx.prior)):
            if ((Nc[i,i] == 0) and (snx.sig0[i] != 0)):
                Nc[i,i] = 1 / snx.sig0[i]**2

        # Remove from the normal matrix of constraints the a priori parameters which do not correspond to any estimated parameter
        # (i.e. RX, RY, RZ in the EMR solutions)
        Nc = Nc[numpy.ix_(ind2, ind2)]

        # Make normal matrix of constraints consistent with the indexing of estimated parameters
        snx.Nc = numpy.zeros((snx.npar, snx.npar))
        snx.Nc[numpy.ix_(ind1,ind1)] = Nc
        Nc = None

        # Total normal matrix
        if (type(snx.Ntot).__name__ != 'NoneType'):
            Ntot = snx.Ntot
            snx.Ntot = None
        elif (type(snx.Q).__name__ != 'NoneType'):
            Ntot = syminv(snx.Q)
            snx.Q = None

        # Unconstrained normal matrix
        snx.N = Ntot - snx.Nc

        # Update prior of sinex object: a priori values are unchanged when available,
        # but set to estimated values when unavailable.
        snx.prior = copy.deepcopy(snx.param)
        snx.sig0  = numpy.zeros(snx.npar)
        for i in range(snx.npar):
            if (snx.Nc[i,i] > 0):
                snx.sig0[i] = 1 / sqrt(snx.Nc[i,i])
        new_x0 = snx.x.copy()
        new_x0[ind1] = snx.x0
        snx.x0 = new_x0

        # Second member of normal equation: k = Ntot * (x-x0)
        snx.k = dot(Ntot, snx.x - snx.x0)

#-------------------------------------------------------------------------------
# Routine : prior2ref
# Purpose : Change a priori values to reference values in a normal equation
# Author  : P. Rebischung
# Created : 11-Aug-2011
#
# Changes :

#
# Input   : - ref       : sinex object containing reference a priori values
#           - log       : Log file. Default is None.
#           - append    : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def prior2ref(self, ref, log=None, append=False):


        # snx = self for class methof # Added JMN 05/01/2018
        snx = self

        # If necessary, open log file and write header
        if (log):
            if (append):
                fl = open(log, 'a')
            else:
                fl = open(log, 'w')

            # Write header
            fl.write('================================================================================\n')
            fl.write('sinex::prior2ref : Change a priori values to reference values in normal equation\n')
            fl.write('================================================================================\n')
            fl.write('\n')

        # Initializations
        dx0 = numpy.zeros(snx.npar)
        b   = False

        # Loop over parameters in normal equation
        for i in range(snx.npar):
            p = snx.param[i]

            # 1 - Case of a SATA parameter
            if (p.type[0:4] == 'SATA'):

                # If this parameter has a reference value,
                if (p.type+p.code in [pref.type+pref.code for pref in ref.param]):
                    iref = [pref.type+pref.code for pref in ref.param].index(p.type+p.code)

                    # And if its reference value is different from its a priori value,
                    if (ref.x[iref] != snx.x0[i]):
                        b = True
                        dx0[i] = ref.x[iref] - snx.x0[i]
                        if (log):
                            fl.write('{0.type} {0.code} {0.pt} {0.soln} : {1:21.14e} => {2:21.14e} ({3:21.14e})\n'.format(p, snx.x0[i], ref.x[iref], ref.x[iref]-snx.x0[i]))

            # 2 - Case of a STA or VEL parameter
            elif (p.type[0:3] in ['STA', 'VEL']):

                # If this parameter has a reference value,
                if (p.type+p.code+p.pt+p.soln in [pref.type+pref.code+pref.pt+pref.soln for pref in ref.param]):
                    iref = [pref.type+pref.code+pref.pt+pref.soln for pref in ref.param].index(p.type+p.code+p.pt+p.soln)

                    # And if its reference value is different from its a priori value,
                    if (ref.x[iref] != snx.x0[i]):
                        b = True
                        dx0[i] = ref.x[iref] - snx.x0[i]
                        if (log):
                            fl.write('{0.type} {0.code} {0.pt} {0.soln} : {1:21.14e} => {2:21.14e} ({3:21.14e})\n'.format(p, snx.x0[i], ref.x[iref], ref.x[iref]-snx.x0[i]))

            # 3 - Case of an ERP
            elif (p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']):

                # If this parameter has a reference value,
                if (p.type+p.tref in [pref.type+pref.tref for pref in ref.param]):
                    iref = [pref.type+pref.tref for pref in ref.param].index(p.type+p.tref)

                    # And if its reference value is different from its a priori value,
                    if (ref.x[iref] != snx.x0[i]):
                        b = True
                        dx0[i] = ref.x[iref] - snx.x0[i]
                        if (log):
                            fl.write('{0.type} {0.code} {0.pt} {0.soln} : {1:21.14e} => {2:21.14e} ({3:21.14e})\n'.format(p, snx.x0[i], ref.x[iref], ref.x[iref]-snx.x0[i]))

            # 4 - Case of a space tie
            elif (p.type in ['TDXJA2', 'TDYJA2', 'TDZJA2', 'TGXJA2', 'TGYJA2', 'TGZJA2', 'TLXJA2', 'TLYJA2', 'TLZJA2']):

                # If this parameter has a reference value,
                if (p.type in [pref.type for pref in ref.param]):
                    iref = [pref.type for pref in ref.param].index(p.type)

                    # And if its reference value is different from its a priori value,
                    if (ref.x[iref] != snx.x0[i]):
                        b = True
                        dx0[i] = ref.x[iref] - snx.x0[i]
                        if (log):
                            fl.write('{0.type} {0.code} {0.pt} {0.soln} : {1:21.14e} => {2:21.14e} ({3:21.14e})\n'.format(p, snx.x0[i], ref.x[iref], ref.x[iref]-snx.x0[i]))


        # If any a priori value has to be modified,
        if (b):

            # Modify a priori values
            snx.x0 = snx.x0 + dx0

            # Modify 2nd member of normal equation
            if (type(snx.k).__name__ != 'NoneType'):
                snx.k = snx.k - dot(snx.N, dx0)

        # Close log file if necessary
        if (log):
            fl.write('\n')
            fl.close()

#-------------------------------------------------------------------------------
# Routine : fix_params
# Purpose : Fix certain types of parameters in a normal equation
# Author  : P. Rebischung
# Created : 11-Aug-2011
#
# Changes :
#
# Input   : types : List of keywords to indicate on which parameters constraints should be removed.
#                   This can include the following keywords:
#                   'XPO', 'XPOR', 'YPO', 'YPOR', 'UT', 'LOD'; 'ERP' for all kinds of ERPs;
#                   'SATA_X', 'SATA_Y', 'SATA_Z'; 'SATA' for all satellite parameters;
#                   'STA'; 'VEL'; 'GC'; 'SBIAS'; 'ST'
# Output  :
#-------------------------------------------------------------------------------
    def fix_params(self, types):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Get indices of parameters to keep
        indk = []
        for i in range(snx.npar):
            p = snx.param[i]
            if not(((p.type == 'XPO   ') and ('XPO'  in types))  or
                   ((p.type == 'XPOR  ') and ('XPOR'  in types)) or
                   ((p.type == 'YPO   ') and ('YPO'  in types))  or
                   ((p.type == 'YPOR  ') and ('YPOR'  in types)) or
                   ((p.type == 'UT    ') and ('UT'  in types))   or
                   ((p.type == 'LOD   ') and ('LOD'  in types))  or
                   ((p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']) and ('ERP' in types)) or
                   ((p.type == 'SATA_X') and ('SATA_X' in types)) or
                   ((p.type == 'SATA_Y') and ('SATA_Y' in types)) or
                   ((p.type == 'SATA_Z') and ('SATA_Z' in types)) or
                   ((p.type[0:4] == 'SATA') and ('SATA' in types)) or
                   ((p.type[0:3] == 'STA') and ('STA' in types)) or
                   ((p.type[0:3] == 'VEL') and ('VEL' in types)) or
                   ((p.type[1:3] == 'GC') and ('GC' in types)) or
                   ((p.type == 'SBIAS ') and ('SBIAS' in types)) or
                   ((p.type in ['TDXJA2', 'TDYJA2', 'TDZJA2', 'TGXJA2', 'TGYJA2', 'TGZJA2', 'TLXJA2', 'TLYJA2', 'TLZJA2']) and ('ST' in types))):
                indk.append(i)

        # Update sinex object
        if (len(indk) < snx.npar):
            snx.npar  = len(indk)
            snx.param = [snx.param[i] for i in indk]
            snx.prior = [snx.prior[i] for i in indk]
            snx.x     = snx.x[indk]
            snx.x0    = snx.x0[indk]
            snx.sig   = snx.sig[indk]
            snx.sig0  = snx.sig0[indk]
            snx.k     = snx.k[indk]
            snx.N     = snx.N[numpy.ix_(indk,indk)]
            if (snx.Nc != None):
                snx.Nc = snx.Nc[numpy.ix_(indk,indk)]

#-------------------------------------------------------------------------------
# Routine : fix_sta
# Purpose : Fix specified stations in a normal equation
# Author  : P. Rebischung
# Created : 11-Aug-2011
#
# Changes :
#
# Input   : - code   : 4-char codes
#           - pt     : PT codes. Default is None.
#           - soln   : Solns. Default is None.
#           - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def fix_sta(self, code, pt=None, soln=None, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # If necessary, open log file and write header
        if (log):
            if (append):
                fl = open(log, 'a')
            else:
                fl = open(log, 'w')
            fl.write('================================================================================\n')
            fl.write('sinex::fix_sta : Fix specified stations from a normal equation\n')
            fl.write('================================================================================\n')
            fl.write('\n')

        # Initialization
        indk = []

        # Loop over parameters
        for i in range(snx.npar):
            p = snx.param[i]

            # Must parameter i be kept?
            keep = False
            if ((pt) and (soln)):
                if (not(p.code+p.pt+p.soln in [code[j]+pt[j]+soln[j] for j in range(len(code))])):
                    keep = True
                    indk.append(i)
            elif (pt):
                if (not(p.code+p.pt in [code[j]+pt[j] for j in range(len(code))])):
                    keep = True
                    indk.append(i)
            else:
                if (not(p.code in code)):
                    keep = True
                    indk.append(i)

            # If parameter i has to be fixed, print message in log file.
            if ((not(keep)) and (log)):
                fl.write('{0.type} {0.code} {0.pt} {0.soln} fixed\n'.format(p))

        # Indices of fixed parameters
        indr = numpy.setdiff1d(list(range(snx.npar)), indk)

        # If any parameter has to be fixed,
        if (len(indr) > 0):

            # Update normal equation
            snx.npar  = len(indk)
            snx.param = [snx.param[i] for i in indk]
            snx.prior = [snx.prior[i] for i in indk]
            snx.x     = snx.x[indk]
            snx.x0    = snx.x0[indk]
            snx.sig   = snx.sig[indk]
            snx.sig0  = snx.sig0[indk]
            snx.N     = snx.N[numpy.ix_(indk,indk)]
            snx.k     = snx.k[indk]
            if (snx.Nc != None):
                snx.Nc = snx.Nc[numpy.ix_(indk,indk)]
            if (snx.Qc != None):
                snx.Qc = snx.Qc[numpy.ix_(indk,indk)]

            # List of points for which there remain estimated positions
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'STAX  ')[0]
            code = [snx.param[i].code for i in ind]
            pt = [snx.param[i].pt for i in ind]
            soln = [snx.param[i].soln for i in ind]

            # Update snx.sta[*].soln
            for s in snx.sta:
                i = 0
                while (i < len(s.soln)):
                    if (not(s.code+s.pt+s.soln[i].soln in [code[k]+pt[k]+soln[k] for k in range(len(code))])):
                        (s.soln).pop(i)
                    else:
                        i = i+1

            # Update snx.sta
            i = 0
            while (i < len(snx.sta)):
                if (len(snx.sta[i].soln) == 0):
                    (snx.sta).pop(i)
                else:
                    i = i+1

        # Close log file if necessary
        if (log):
            fl.write('\n')
            fl.close()

#-------------------------------------------------------------------------------
# Routine : reduce_params
# Purpose : Reduce certain types of parameters from a normal equation
# Author  : P. Rebischung
# Created : 26-Jun-2011
#
# Changes :
#
# Input   : types : List of keywords to indicate on which parameters constraints should be removed.
#                   This can include the following keywords:
#                   'XPO', 'XPOR', 'YPO', 'YPOR', 'UT', 'LOD'; 'ERP' for all kinds of ERPs;
#                   'SATA_X', 'SATA_Y', 'SATA_Z'; 'SATA' for all satellite parameters;
#                   'STA'; 'VEL'; 'GC'; 'RBIAS'
# Output  :
#-------------------------------------------------------------------------------
    def reduce_params(self, types):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Return if no parameter has to be reduced
        if (len(types) == 0):
            return

        # Get indices of parameters to keep
        indk = []
        for i in range(snx.npar):
            p = snx.param[i]
            if not(((p.type == 'XPO   ') and ('XPO'  in types))  or
                   ((p.type == 'XPOR  ') and ('XPOR'  in types)) or
                   ((p.type == 'YPO   ') and ('YPO'  in types))  or
                   ((p.type == 'YPOR  ') and ('YPOR'  in types)) or
                   ((p.type == 'UT    ') and ('UT'  in types))   or
                   ((p.type == 'LOD   ') and ('LOD'  in types))  or
                   ((p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']) and ('ERP' in types)) or
                   ((p.type == 'SATA_X') and ('SATA_X' in types)) or
                   ((p.type == 'SATA_Y') and ('SATA_Y' in types)) or
                   ((p.type == 'SATA_Z') and ('SATA_Z' in types)) or
                   ((p.type[0:4] == 'SATA') and ('SATA' in types)) or
                   ((p.type[0:3] == 'STA') and ('STA' in types)) or
                   ((p.type[0:3] == 'VEL') and ('VEL' in types)) or
                   ((p.type in ['XGC   ', 'YGC   ', 'ZGC   ', 'SC    ']) and ('GC' in types)) or
                   ((p.type == 'RBIAS ') and ('RBIAS' in types))):
                indk.append(i)

        # Indices of reduced parameters
        indr = numpy.setdiff1d(list(range(snx.npar)), indk)

        # Update sinex object
        if (len(indr) > 0):
            snx.npar  = len(indk)
            snx.param = [snx.param[i] for i in indk]
            snx.prior = [snx.prior[i] for i in indk]
            snx.x     = snx.x[indk]
            snx.x0    = snx.x0[indk]
            snx.sig   = snx.sig[indk]
            snx.sig0  = snx.sig0[indk]
            R         = dot(snx.N[numpy.ix_(indk,indr)], syminv(snx.N[numpy.ix_(indr,indr)]))
            snx.N     = snx.N[numpy.ix_(indk,indk)] - dot(R, snx.N[numpy.ix_(indr,indk)])
            snx.k     = snx.k[indk] - dot(R, snx.k[indr])
            if (snx.Nc != None):
                snx.Nc = snx.Nc[numpy.ix_(indk,indk)]

#-------------------------------------------------------------------------------
# Routine : reduce_sta
# Purpose : Reduce specified stations from a normal equation
# Author  : P. Rebischung
# Created : 02-Sep-2014
#
# Changes :
#
# Input   : - code   : 4-char codes
#           - pt     : PT codes. Default is None.
#           - soln   : Solns. Default is None.
#           - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def reduce_sta(self, code, pt=None, soln=None, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # If necessary, open log file and write header
        if (log):
            if (append):
                fl = open(log, 'a')
            else:
                fl = open(log, 'w')
            fl.write('================================================================================\n')
            fl.write('sinex::reduce_sta : Reduce specified stations from a normal equation\n')
            fl.write('================================================================================\n')
            fl.write('\n')

        # Initialization
        indk = []

        # Loop over parameters
        for i in range(snx.npar):
            p = snx.param[i]

            # Must parameter i be kept?
            keep = False
            if ((pt) and (soln)):
                if (not(p.code+p.pt+p.soln in [code[j]+pt[j]+soln[j] for j in range(len(code))])):
                    keep = True
                    indk.append(i)
            elif (pt):
                if (not(p.code+p.pt in [code[j]+pt[j] for j in range(len(code))])):
                    keep = True
                    indk.append(i)
            else:
                if (not(p.code in code)):
                    keep = True
                    indk.append(i)

            # If parameter i has to be removed, print message in log file.
            if ((not(keep)) and (log)):
                fl.write('{0.type} {0.code} {0.pt} {0.soln} reduced\n'.format(p))

        # Indices of reduced parameters
        indr = numpy.setdiff1d(list(range(snx.npar)), indk)

        # If any parameter has to be reduced,
        if (len(indr) > 0):

            # Update normal equation
            snx.npar  = len(indk)
            snx.param = [snx.param[i] for i in indk]
            snx.prior = [snx.prior[i] for i in indk]
            snx.x     = snx.x[indk]
            snx.x0    = snx.x0[indk]
            snx.sig   = snx.sig[indk]
            snx.sig0  = snx.sig0[indk]
            R         = dot(snx.N[numpy.ix_(indk,indr)], syminv(snx.N[numpy.ix_(indr,indr)]))
            snx.N     = snx.N[numpy.ix_(indk,indk)] - dot(R, snx.N[numpy.ix_(indr,indk)])
            snx.k     = snx.k[indk] - dot(R, snx.k[indr])
            if (snx.Nc != None):
                snx.Nc = snx.Nc[numpy.ix_(indk,indk)]

            # List of points for which there remain estimated positions
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'STAX  ')[0]
            code = [snx.param[i].code for i in ind]
            pt = [snx.param[i].pt for i in ind]
            soln = [snx.param[i].soln for i in ind]

            # Update snx.sta[*].soln
            for s in snx.sta:
                i = 0
                while (i < len(s.soln)):
                    if (not(s.code+s.pt+s.soln[i].soln in [code[k]+pt[k]+soln[k] for k in range(len(code))])):
                        (s.soln).pop(i)
                    else:
                        i = i+1

            # Update snx.sta
            i = 0
            while (i < len(snx.sta)):
                if (len(snx.sta[i].soln) == 0):
                    (snx.sta).pop(i)
                else:
                    i = i+1

        # Close log file if necessary
        if (log):
            fl.write('\n')
            fl.close()

#-------------------------------------------------------------------------------
# Routine : reduce_doublons
# Purpose : Reduce stations appearing twice in a normal equation
# Author  : P. Rebischung
# Created : 21-Aug-2013
#
# Changes :
#
# Input   : - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def reduce_doublons(self, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initializations
        code = []
        pt   = []

        # Loop over STAX parameters
        for p in snx.param:
            if (p.type == 'STAX  '):

                # If current stations appears twice, add it to the red list if it is not yet in.
                ind = numpy.nonzero(numpy.array([q.type+q.code+q.pt for q in snx.param]) == 'STAX  '+p.code+p.pt)[0]
                if ((len(ind) > 1) and (not(p.code+p.pt in [code[i]+pt[i] for i in range(len(code))]))):
                    code.append(p.code)
                    pt.append(p.pt)

        # Reduce stations appearing twice
        snx.reduce_sta(code, pt, log=log, append=append)

#-------------------------------------------------------------------------------
# Routine : setup_geocenter
# Purpose : Add geocenter coordinates into a normal equation
# Author  : P. Rebischung
# Created : 02-Sep-2014
#
# Changes :
#
# Input   : tref : Reference epoch (SINEX format)
# Output  :
#-------------------------------------------------------------------------------
    def setup_geocenter(self, tref):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Add XGC parameter
        p = record()
        p.type  = 'XGC   '
        p.code  = '----'
        p.pt    = '--'
        p.soln  = '----'
        p.tref  = tref
        p.unit  = 'm   '
        p.const = '2'
        snx.param.append(p)
        snx.prior.append(p)

        # Add YGC parameter
        p = copy.deepcopy(p)
        p.type = 'YGC   '
        snx.param.append(p)
        snx.prior.append(p)

        # Add ZGC parameter
        p = copy.deepcopy(p)
        p.type = 'ZGC   '
        snx.param.append(p)
        snx.prior.append(p)

        # Add a priori geocenter coordinates
        snx.x0 = numpy.hstack((snx.x0, [0, 0, 0]))
        snx.sig0 = numpy.hstack((snx.sig0, [0, 0, 0]))

        # Update snx.Nc
        if (snx.Nc != None):
            snx.Nc = numpy.vstack((numpy.hstack((snx.Nc, numpy.zeros((snx.npar, 3)))), numpy.zeros((3, snx.npar+3))))

        # "Design" matrix
        A = numpy.zeros((snx.npar+3, snx.npar))
        A[0:snx.npar, 0:snx.npar] = numpy.eye(snx.npar)
        for i in range(snx.npar):
            if (snx.param[i].type == 'STAX  '):
                A[snx.npar:snx.npar+3, i:i+3] = -numpy.eye(3)

        # Update snx.N and snx.k
        snx.N = dot(A, dot(snx.N, A.T))
        snx.k = dot(A, snx.k)

        # Update snx.npar
        snx.npar = snx.npar + 3

#-------------------------------------------------------------------------------
# Routine : setup_scale
# Purpose : Add terrestrial scale into a normal equation
# Author  : P. Rebischung
# Created : 02-Sep-2014
#
# Changes :
#
# Input   : tref : Reference epoch (SINEX format)
# Output  :
#-------------------------------------------------------------------------------
    def setup_scale(self, tref):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Add SC parameter
        p = record()
        p.type  = 'SC    '
        p.code  = '----'
        p.pt    = '--'
        p.soln  = '----'
        p.tref  = tref
        p.unit  = 'ppb '
        p.const = '2'
        snx.param.append(p)
        snx.prior.append(p)

        # Add a priori scale factor
        snx.x0 = numpy.hstack((snx.x0, [0]))
        snx.sig0 = numpy.hstack((snx.sig0, [0]))

        # Update snx.Nc
        if (snx.Nc != None):
            snx.Nc = numpy.vstack((numpy.hstack((snx.Nc, numpy.zeros((snx.npar, 1)))), numpy.zeros((1, snx.npar+1))))

        # "Design" matrix
        A = numpy.zeros((snx.npar+1, snx.npar))
        A[0:snx.npar, 0:snx.npar] = numpy.eye(snx.npar)
        for i in range(snx.npar):
            if (snx.param[i].type[0:3] == 'STA'):
                A[snx.npar, i] = -1e-9*snx.x[i]

        # Update snx.N and snx.k
        snx.N = dot(A, dot(snx.N, A.T))
        snx.k = dot(A, snx.k)

        # Update snx.npar
        snx.npar = snx.npar + 1

#-------------------------------------------------------------------------------
# Routine : reduce_helmerts
# Purpose : Reduce origin, scale and/or orientation in a normal equation
# Author  : P. Rebischung
# Created : 31-Jul-2013
#
# Changes :
#
# Input   : params : String indicating which Helmert parameters should be reduced.
#                    It can be composed by any combination of letters 'T' (translations),
#                   'S' (scale) and 'R' (rotations).
# Output  :
#-------------------------------------------------------------------------------
    def reduce_helmerts(self, params):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Get Helmert design matrix
        A = snx.get_helmert_matrix(params, pole=True)

        # Non-orthogonal projector onto Ker(A')
        ANAi = sympinv(dot(A.T, dot(snx.N, A)))
        P    = numpy.eye(snx.npar) - dot(snx.N, dot(A, dot(ANAi, A.T)))

        # Update normal equation
        snx.N = dot(P, snx.N)
        snx.k = dot(P, snx.k)

#-------------------------------------------------------------------------------
# Routine : add_min_const
# Purpose : Add minimal constraints to normal matrix of constraints (snx.Nc)
# Author  : P. Rebischung
# Created : 08-Mar-2012
#
# Changes :
#
# Input   : - params         : String indicating on which parameters minimal constraints should be applied.
#                              It can be composed by any combination of letters 'T' (translations),
#                              'S' (scale) and 'R' (rotations). Default is '' (no minimal constraints).
#           - sigma          : Sigma of minimal constraints in m. Default is 1e-4.
#           - code, pt, soln : If specified, then minimal constraints will be applied to these points only.
#           - reduce_B       : True if the rows of B should be reduced rather than the columns of A. Default is False.
# Output  :
#-------------------------------------------------------------------------------
    def add_min_const(self, params='', sigma=1e-4, code=None, pt=None, soln=None, reduce_B=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Loop over STAX parameters in order to compute design matrix A of minimal constraints
        A = numpy.zeros((snx.npar, 7))
        for i in range(snx.npar):
            p = snx.param[i]
            if (p.type == 'STAX  '):

                # Check whether minimal constraints have to be applied to current point
                b = False
                if ((code) and (pt) and (soln)):
                    if (p.code+p.pt+p.soln in [code[j]+pt[j]+soln[j] for j in range(len(code))]):
                        b = True
                elif ((code) and (pt)):
                    if (p.code+p.pt in [code[j]+pt[j] for j in range(len(code))]):
                        b = True
                elif (code):
                    if (p.code in [code[j] for j in range(len(code))]):
                        b = True
                else:
                    b = True

                # If minimal constraints have to be applied to current point, fill design matrix
                if (b):
                    A[i+0, 0] =  1
                    A[i+1, 1] =  1
                    A[i+2, 2] =  1
                    A[i+0, 3] =  snx.x0[i+0] / ae
                    A[i+1, 3] =  snx.x0[i+1] / ae
                    A[i+2, 3] =  snx.x0[i+2] / ae
                    A[i+1, 4] = -snx.x0[i+2] / ae
                    A[i+2, 4] =  snx.x0[i+1] / ae
                    A[i+0, 5] =  snx.x0[i+2] / ae
                    A[i+2, 5] = -snx.x0[i+0] / ae
                    A[i+0, 6] = -snx.x0[i+1] / ae
                    A[i+1, 6] =  snx.x0[i+0] / ae

        # Indices of Helmert parameters that should effectively be constrained
        ind = []
        if ('T' in params):
            ind.extend(list(range(0, 3)))
        if ('S' in params):
            ind.append(3)
        if ('R' in params):
            ind.extend(list(range(4, 7)))

        # 1st case: reduce columns of A and compute B
        if (not(reduce_B)):
            A = A[:,ind]
            B = dot(syminv(dot(A.T, A)), A.T)

        # 2nd case: compute B and reduce rows of B
        else:
            B = dot(syminv(dot(A.T, A)), A.T)
            B = B[ind,:]

        # Add minimal constraints to normal matrix of constraints
        snx.Nc = snx.Nc + dot(B.T, B) / sigma**2

#-------------------------------------------------------------------------------
# Routine : set_constraints
# Purpose : Define normal matrix of constraints (snx.Nc)
# Author  : P. Rebischung
# Created : 08-Mar-2012
#
# Changes :
#
# Input   : - keep_const_on  : List of keywords indicating for which parameters constraints should be kept.
#                              This can include the following keywords: 'STA' for all station positions; 'VEL' for all station velocities;
#                              'XPO', 'XPOR', 'YPO', 'YPOR', 'UT', 'LOD'; 'ERP' for all kinds of ERPs;
#                              'SATA_X', 'SATA_Y', 'SATA_Z'; 'SATA' for all satellite antenna parameters and 'GC'.
#                              Default is [] (no constraints kept).
#           - add_const_on   : List of keywords indicating for which parameters constraints should be added.
#                              It can include the same keywords as keep_const_on except 'ERP'.
#                              Default is [] (no constraints added).
#           - add_const_sig  : Sigmas of constraints which should be added. (One value per keyword in add_const_on.)
#                              Default is [].
#           - min_const_on   : String indicating on which parameters minimal constraints should be applied.
#                              It can be composed by any combination of letters 'T' (translations),
#                              'S' (scale) and 'R' (rotations). Default is '' (no minimal constraints).
#           - min_const_sig  : Sigma of minimal constraints [m].
#           - code, pt, soln : If specified, then minimal constraints will be applied to these points only.
#           - reduce_B       : True if the rows of B should be reduced rather than the columns of A. Default is False.
# Output  :
#-------------------------------------------------------------------------------
    def set_constraints(self, keep_const_on=[], add_const_on=[], add_const_sig=[], min_const_on='',
                        min_const_sig=1e-4, code=None, pt=None, soln=None, reduce_B=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # 1 - Re-initialize normal matrix of constraints (snx.Nc)
        #--------------------------------------------------------

        # If none of the original constraints should be kept, then set snx.Nc = 0
        if (len(keep_const_on) == 0):
            snx.Nc = numpy.zeros((len(snx.param), len(snx.param)))

        # Else,
        else:

            # First, if snx.Nc is not defined, define it as a diagonal matrix
            # from the standard deviations of the a priori parameters.
            if (snx.Nc is None):
                sc      = snx.sig0
                ind     = numpy.nonzero(sc == 0)[0]
                sc[ind] = numpy.inf
                snx.Nc  = numpy.diag(1 / sc**2)

            # Get indices of parameters for which the original constraints should be kept
            indc = []
            for i in range(len(snx.param)):
                p = snx.param[i]
                if (((p.type[0:3] == 'STA') and ('STA' in keep_const_on)) or
                    ((p.type[0:3] == 'VEL') and ('VEL' in keep_const_on)) or
                    ((p.type == 'XPO   ') and ('XPO'  in keep_const_on))  or
                    ((p.type == 'XPOR  ') and ('XPOR'  in keep_const_on)) or
                    ((p.type == 'YPO   ') and ('YPO'  in keep_const_on))  or
                    ((p.type == 'YPOR  ') and ('YPOR'  in keep_const_on)) or
                    ((p.type == 'UT    ') and ('UT'  in keep_const_on))   or
                    ((p.type == 'LOD   ') and ('LOD'  in keep_const_on))  or
                    ((p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']) and ('ERP' in keep_const_on)) or
                    ((p.type == 'SATA_X') and ('SATA_X' in keep_const_on)) or
                    ((p.type == 'SATA_Y') and ('SATA_Y' in keep_const_on)) or
                    ((p.type == 'SATA_Z') and ('SATA_Z' in keep_const_on)) or
                    ((p.type[0:4] == 'SATA') and ('SATA' in keep_const_on)) or
                    ((p.type[1:3] == 'GC') and ('GC' in keep_const_on))):

                    indc.append(i)

            # Keep snx.Nc[indc,indc] unchanged and set all other elements to 0
            Nc = snx.Nc[numpy.ix_(indc, indc)]
            snx.Nc = numpy.zeros((len(snx.param), len(snx.param)))
            snx.Nc[numpy.ix_(indc, indc)] = snx.Nc[numpy.ix_(indc, indc)] + Nc

        # Update standard deviations of a priori parameters
        snx.sig0 = numpy.zeros(snx.npar)
        ind = numpy.nonzero(numpy.diag(snx.Nc) != 0)
        snx.sig0[ind] = 1 / numpy.sqrt(numpy.diag(snx.Nc)[ind])

        # 2 - Add diagonal constaints if requested
        #-----------------------------------------

        # On station coordinates
        if ('STA' in add_const_on):
            ic = add_const_on.index('STA')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type[0:3] for p in snx.param]) == 'STA')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On XPO
        if ('XPO' in add_const_on):
            ic = add_const_on.index('XPO')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'XPO   ')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On XPOR
        if ('XPOR' in add_const_on):
            ic = add_const_on.index('XPOR')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'XPOR  ')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On YPO
        if ('YPO' in add_const_on):
            ic = add_const_on.index('YPO')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'YPO   ')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On YPOR
        if ('YPOR' in add_const_on):
            ic = add_const_on.index('YPOR')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'YPOR  ')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On UT
        if ('UT' in add_const_on):
            ic = add_const_on.index('UT')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'UT    ')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On LOD
        if ('LOD' in add_const_on):
            ic = add_const_on.index('LOD')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'LOD   ')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On SATA_X
        if ('SATA_X' in add_const_on):
            ic = add_const_on.index('SATA_X')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'SATA_X')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On SATA_Y
        if ('SATA_Y' in add_const_on):
            ic = add_const_on.index('SATA_Y')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'SATA_Y')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On SATA_Z
        if ('SATA_Z' in add_const_on):
            ic = add_const_on.index('SATA_Z')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'SATA_Z')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On all SATA parameters
        if ('SATA' in add_const_on):
            ic = add_const_on.index('SATA')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type[0:4] for p in snx.param]) == 'SATA')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # On geocenter
        if ('GC' in add_const_on):
            ic = add_const_on.index('GC')
            wc = 1 / add_const_sig[ic]**2
            ind = numpy.nonzero(numpy.array([p.type[1:3] for p in snx.param]) == 'GC')
            snx.Nc[ind, ind] = snx.Nc[ind, ind] + wc

        # 3 - Add minimal constraints if requested
        #-----------------------------------------
        if (len(min_const_on) > 0):
            snx.add_min_const(params=min_const_on, sigma=min_const_sig, code=code, pt=pt, soln=soln, reduce_B=reduce_B)

#-------------------------------------------------------------------------------
# Routine : neqinv
# Purpose : Invert normal equation
# Author  : P. Rebischung
# Created : 07-Jul-2011
#
# Changes :
#
# Input   : - keep_const_on  : List of keywords indicating for which parameters constraints should be kept.
#                              This can include the following keywords: 'STA' for all station positions; 'VEL' for all station velocities;
#                              'XPO', 'XPOR', 'YPO', 'YPOR', 'UT', 'LOD'; 'ERP' for all kinds of ERPs;
#                              'SATA_X', 'SATA_Y', 'SATA_Z'; 'SATA' for all satellite antenna parameters and 'GC'.
#                              Default is [] (no constraints kept).
#           - add_const_on   : List of keywords indicating for which parameters constraints should be added.
#                              It can include the same keywords as keep_const_on except 'ERP'.
#                              Default is [] (no constraints added).
#           - add_const_sig  : Sigmas of constraints which should be added. (One value per keyword in add_const_on.)
#                              Default is [].
#           - min_const_on   : String indicating on which parameters minimal constraints should be applied.
#                              It can be composed by any combination of letters 'T' (translations),
#                              'S' (scale) and 'R' (rotations). Default is '' (no minimal constraints).
#           - min_const_sig  : Sigma of minimal constraints [m].
#           - code, pt, soln : If specified, then minimal constraints will be applied to these points only.
#           - reduce_B       : True if the rows of B should be reduced rather than the columns of A. Default is False.
# Output  :
#-------------------------------------------------------------------------------
    def neqinv(self, keep_const_on=[], add_const_on=[], add_const_sig=[], min_const_on='',
               min_const_sig=1e-4, code=None, pt=None, soln=None, reduce_B=False):

        
        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Set normal matrix of constraints
        snx.set_constraints(keep_const_on, add_const_on, add_const_sig, min_const_on, min_const_sig, code, pt, soln, reduce_B)

        # Solve normal equation
        snx.Q = syminv(snx.N + snx.Nc)
        snx.x = snx.x0 + dot(snx.Q, snx.k)
        snx.sig = numpy.sqrt(numpy.diag(snx.Q))

        # Clear snx.N and snx.k
        snx.N = None
        snx.k = None

#-------------------------------------------------------------------------------
# Routine : rescale
# Purpose : Multiply solution covariance matrix by given factor
# Author  : P. Rebischung
# Created : 20-May-2012
#
# Changes :
#
# Input   : f : Scale factor (square-root of variance factor)
# Output  :
#-------------------------------------------------------------------------------
    def rescale(self, f):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Rescale sigmas
        snx.sig = snx.sig * f

        # Rescals covariance matrix
        if (snx.Q != None):
            snx.Q = f**2 * snx.Q

#-------------------------------------------------------------------------------
# Routine : get_common_points
# Purpose : Get list of common points between two solutions
# Author  : P. Rebischung
# Created : 09-Mar-2012
#
# Changes :
#
# Input   : - ref       : Second sinex object
# Output  : - code      : 4-char codes of common points
#           - pt        : PT codes of common points
#           - soln      : Solns of common points
#-------------------------------------------------------------------------------
    def get_common_points(self, ref):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initializations
        code = []
        pt = []
        soln = []

        # Loop over STAX parameters of snx to identify common stations
        for i in range(snx.npar):
            p1 = snx.param[i]
            if (p1.type == 'STAX  '):

                # If current station is common to both solutions,
                if ('STAX  '+p1.code+p1.pt+p1.soln in [p0.type+p0.code+p0.pt+p0.soln for p0 in ref.param]):
                    code.append(p1.code)
                    pt.append(p1.pt)
                    soln.append(p1.soln)

        return (code, pt, soln)

#-------------------------------------------------------------------------------
# Routine : helmert_wrt
# Purpose : Helmert comparison between two solutions
# Author  : P. Rebischung
# Created : 24-May-2011
#
# Changes :
#
# Input   : - ref       : Solution to which the comparison is made (sinex object)
#           - params    : String indicating which parameters should be estimated.
#                         It can be composed by any combination of letters 'T' (translations),
#                         'S' (scale) and 'R' (rotations). Default is 'RST' (all 7 parameters).
#           - weighting : Keyword to indicate which weighting should be used.
#                         It can take the following values :
#                         - 'identity' to use an identity weight matrix (default)
#                         - 'diagonal' to use a diagonal weight matrix
#                         - 'full' to use a full weight matrix (inv(snx.Q))
#           - with_vel  : True or False depending on whether 7 or 14 parameters
#                         must be estimated. Default is False.
#           - normalize : If False, covariance matrix of estimated Helmert parameters
#                         will not be scaled by unit variance factor. Default is True.
#           - log       : Log file. Default is None.
#           - append    : If true, text will be appended to log file. Default if False.
# Output  : - T         : 7 (14) Helmert parameters (mm, ppb, mas, [mm/y, ppb/y, mas/y])
#           - Q         : Corresponding covariance matrix
#           - code      : 4-char codes of points used for the comparison
#           - pt        : PT codes of points used for the comparison
#           - soln      : solns of points used for the comparison
#           - vl        : E, N, H residuals (mm[/y])
#           - vn        : Normalized E, N, H residuals
#           - w         : E, N, H WRMS (mm[/y])
#-------------------------------------------------------------------------------
    def helmert_wrt(self, ref, params='RST', weighting='identity', with_vel=False, normalize=True, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initializations
        nsta = 0
        ind1 = []
        ind0 = []
        code = []
        pt   = []
        soln = []

        # Loop over STAX parameters of snx to identify common stations
        for i in range(len(snx.param)):
            p1 = snx.param[i]
            if (p1.type == 'STAX  '):

                # If current station is common to both solutions,
                if ('STAX  '+p1.code+p1.pt+p1.soln in [p0.type+p0.code+p0.pt+p0.soln for p0 in ref.param]):
                    nsta = nsta+1
                    code.append(p1.code)
                    pt.append(p1.pt)
                    soln.append(p1.soln)

                    j = [p0.type+p0.code+p0.pt+p0.soln for p0 in ref.param].index('STAX  '+p1.code+p1.pt+p1.soln)
                    if (with_vel):
                        ind1.extend(list(range(i, i+6)))
                        ind0.extend(list(range(j, j+6)))
                    else:
                        ind1.extend(list(range(i, i+3)))
                        ind0.extend(list(range(j, j+3)))

        # Get reference coordinates of common stations
        nobs = len(ind0)
        X = numpy.zeros((nsta, 3))
        if (with_vel):
            X[:,0] = ref.x[ind0[0:nobs:6]]
            X[:,1] = ref.x[ind0[1:nobs:6]]
            X[:,2] = ref.x[ind0[2:nobs:6]]
        else:
            X[:,0] = ref.x[ind0[0:nobs:3]]
            X[:,1] = ref.x[ind0[1:nobs:3]]
            X[:,2] = ref.x[ind0[2:nobs:3]]

        # Design matrix in case of 14-parameter comparison
        if (with_vel):
            A = numpy.zeros((nobs, 14))

            A[0:nobs:6, 0]  =  1e-3
            A[1:nobs:6, 1]  =  1e-3
            A[2:nobs:6, 2]  =  1e-3
            A[0:nobs:6, 3]  =  1e-9 * X[:,0]
            A[1:nobs:6, 3]  =  1e-9 * X[:,1]
            A[2:nobs:6, 3]  =  1e-9 * X[:,2]
            A[1:nobs:6, 4]  = -mas2rad * X[:,2]
            A[2:nobs:6, 4]  =  mas2rad * X[:,1]
            A[0:nobs:6, 5]  =  mas2rad * X[:,2]
            A[2:nobs:6, 5]  = -mas2rad * X[:,0]
            A[0:nobs:6, 6]  = -mas2rad * X[:,1]
            A[1:nobs:6, 6]  =  mas2rad * X[:,0]

            A[3:nobs:6, 7]  =  1e-3
            A[4:nobs:6, 8]  =  1e-3
            A[5:nobs:6, 9]  =  1e-3
            A[3:nobs:6, 10] =  1e-9 * X[:,0]
            A[4:nobs:6, 10] =  1e-9 * X[:,1]
            A[5:nobs:6, 10] =  1e-9 * X[:,2]
            A[4:nobs:6, 11] = -mas2rad * X[:,2]
            A[5:nobs:6, 11] =  mas2rad * X[:,1]
            A[3:nobs:6, 12] =  mas2rad * X[:,2]
            A[5:nobs:6, 12] = -mas2rad * X[:,0]
            A[3:nobs:6, 13] = -mas2rad * X[:,1]
            A[4:nobs:6, 13] =  mas2rad * X[:,0]

        # Design matrix in case of 7-parameter comparison
        else:
            A = numpy.zeros((nobs, 7))

            A[0:nobs:3, 0]  =  1e-3
            A[1:nobs:3, 1]  =  1e-3
            A[2:nobs:3, 2]  =  1e-3
            A[0:nobs:3, 3]  =  1e-9 * X[:, 0]
            A[1:nobs:3, 3]  =  1e-9 * X[:, 1]
            A[2:nobs:3, 3]  =  1e-9 * X[:, 2]
            A[1:nobs:3, 4]  = -mas2rad * X[:, 2]
            A[2:nobs:3, 4]  =  mas2rad * X[:, 1]
            A[0:nobs:3, 5]  =  mas2rad * X[:, 2]
            A[2:nobs:3, 5]  = -mas2rad * X[:, 0]
            A[0:nobs:3, 6]  = -mas2rad * X[:, 1]
            A[1:nobs:3, 6]  =  mas2rad * X[:, 0]

        # Indices of parameters that have to be estimated (ind)
        ind  = []

        i = 0
        if ('T' in params):
            ind.extend(list(range(0,3)))
        if ('S' in params):
            ind.append(3)
        if ('R' in params):
            ind.extend(list(range(4,7)))

        if (with_vel):
            if ('T' in params):
                ind.extend(list(range(7,10)))
            if ('S' in params):
                ind.append(10)
            if ('R' in params):
                ind.extend(list(range(11,14)))

        # Reduce A to these parameters
        A = A[:,ind]
        npar = len(ind)

        # 2nd member B
        B = snx.x[ind1] - ref.x[ind0]

        # Compute normal matrix and 2nd member of normal equation
        if (weighting == 'identity'):
            N  = dot(A.T, A)
            K  = dot(A.T, B)
        elif (weighting == 'diagonal'):
            P  = numpy.diag(1 / snx.sig[ind1]**2)
            N  = A.T * P
            K  = dot(N, B)
            N  = dot(N, A)
        elif (weighting == 'full'):
            P  = syminv(snx.Q[numpy.ix_(ind1,ind1)])
            N  = dot(A.T, P)
            K  = dot(N, B)
            N  = dot(N, A)

        # Solve normal equation
        Q = syminv(N)
        T = dot(Q, K)

        # Compute residuals and variance factor
        v = B - dot(A, T)

        if (weighting == 'identity'):
            vPv = dot(v, v)
        elif (weighting == 'diagonal'):
            vPv = dot(v*P, v)
        elif (weighting == 'full'):
            vPv = dot(v, dot(P, v))

        sigma0 = sqrt(vPv / (nobs - npar))

        # Compute covariance matrix of residuals
        if (weighting == 'identity'):
            Qv = numpy.eye(nobs) - dot(A, dot(Q, A.T))
        elif (weighting == 'diagonal'):
            Qv = numpy.diag(snx.sig[ind1]**2) - dot(A, dot(Q, A.T))
        else:
            Qv = snx.Q[numpy.ix_(ind1,ind1)] - dot(A, dot(Q, A.T))

        # Scale covariance matrices with unit variance factor
        if (normalize):
            Q = sigma0**2 * Q
            Qv = sigma0**2 * Qv

        # Compute (X,Y,Z) -> (E,N,H) rotation matrices
        R = xyz2enh(X)

        # Compute E,N,H [normalized] residuals and E,N,H WRMS
        vl = numpy.zeros(nobs)
        vn = numpy.zeros(nobs)

        # 1st case : 14 parameters
        if (with_vel):
            w  = numpy.zeros(6)
            d  = numpy.zeros(6)

            for i in range(nsta):
                vl[6*i+0:6*i+3] = dot(R[i], v[6*i+0:6*i+3])
                Ql = dot(R[i], dot(Qv[6*i+0:6*i+3, 6*i+0:6*i+3], R[i].T))
                vn[6*i+0:6*i+3] = vl[6*i+0:6*i+3] / numpy.sqrt(numpy.diag(Ql))
                w[0:3] = w[0:3] + vl[6*i+0:6*i+3]**2 / numpy.diag(Ql)
                d[0:3] = d[0:3] + 1 / numpy.diag(Ql)

                vl[6*i+3:6*i+6] = dot(R[i], v[6*i+3:6*i+6])
                Ql = dot(R[i], dot(Qv[6*i+3:6*i+6, 6*i+3:6*i+6], R[i].T))
                vn[6*i+3:6*i+6] = vl[6*i+3:6*i+6] / numpy.sqrt(numpy.diag(Ql))
                w[3:6] = w[3:6] + vl[6*i+3:6*i+6]**2 / numpy.diag(Ql)
                d[3:6] = d[3:6] + 1 / numpy.diag(Ql)

        # 2nd case : 7 parameters
        else:
            w  = numpy.zeros(3)
            d  = numpy.zeros(3)

            for i in range(nsta):
                vl[3*i+0:3*i+3] = dot(R[i], v[3*i+0:3*i+3])
                Ql = dot(R[i], dot(Qv[3*i+0:3*i+3, 3*i+0:3*i+3], R[i].T))
                vn[3*i+0:3*i+3] = vl[3*i+0:3*i+3] / numpy.sqrt(numpy.diag(Ql))
                w = w + vl[3*i+0:3*i+3]**2 / numpy.diag(Ql)
                d = d + 1 / numpy.diag(Ql)

        w = numpy.sqrt(w / d)

        # Reshape residual arrays
        if (with_vel):
            vl.resize(nsta, 6)
            vn.resize(nsta, 6)
        else:
            vl.resize(nsta, 3)
            vn.resize(nsta, 3)

        # Convert raw residuals and WRMS into mm[/y]
        vl = 1000 * vl
        w  = 1000 * w

        # Reshape array of transformation parameters and covariance matrix
        # (Put zeros for non-estimated parameters)
        if (with_vel):
            Tfull  = numpy.zeros(14)
            Qfull = numpy.zeros((14, 14))
        else:
            Tfull  = numpy.zeros(7)
            Qfull = numpy.zeros((7, 7))

        Tfull[ind] = T
        Qfull[numpy.ix_(ind,ind)] = Q

        T = Tfull
        Q = Qfull
        sT = numpy.sqrt(numpy.diag(Q))

        # Write log file
        if (log):
            if (append):
                f = open(log, 'a')
            else:
                f = open(log, 'w')

            # Write header
            f.write('================================================================================\n')
            f.write('sinex::helmert_wrt : {0} vs {1}\n'.format(snx.file, ref.file))
            f.write('================================================================================\n')
            f.write('\n')

            # Write main options and statistics
            f.write('Number of common points : {0}\n'.format(nsta))
            f.write('Number of parameters    : {0}\n'.format(npar))
            f.write('Weighting               : {0}\n'.format(weighting))
            f.write('Variance factor         : {0:.8e}\n'.format(sigma0))
            f.write('WRMS East               : {0:8.3f} mm\n'.format(w[0]))
            f.write('WRMS North              : {0:8.3f} mm\n'.format(w[1]))
            f.write('WRMS Up                 : {0:8.3f} mm\n'.format(w[2]))
            if (with_vel):
                f.write('WRMS vel East           : {0:8.3f} mm/y\n'.format(w[3]))
                f.write('WRMS vel North          : {0:8.3f} mm/y\n'.format(w[4]))
                f.write('WRMS vel Up             : {0:8.3f} mm/y\n'.format(w[5]))
            f.write('\n')

            # Write estimated parameters and formal errors
            f.write('Estimated Helmert parameters :\n')
            f.write('------------------------------\n')
            f.write('\n')
            f.write('       |  TX(mm)  TY(mm)  TZ(mm) SC(ppb) RX(mas) RY(mas) RZ(mas) |\n')
            f.write('-------|---------------------------------------------------------|\n')
            f.write('       | {0[0]:7.1f} {0[1]:7.1f} {0[2]:7.1f} {0[3]:7.2f} {0[4]:7.3f} {0[5]:7.3f} {0[6]:7.3f} |\n'.format(T))
            f.write('  +/-  | {0[0]:7.1f} {0[1]:7.1f} {0[2]:7.1f} {0[3]:7.2f} {0[4]:7.3f} {0[5]:7.3f} {0[6]:7.3f} |\n'.format(sT))
            f.write('-------|---------------------------------------------------------|\n')

            if (with_vel):
                f.write(' rates |  (mm/y)  (mm/y)  (mm/y) (ppb/y) (mas/y) (mas/y) (mas/y) |\n')
                f.write('-------|---------------------------------------------------------|\n')
                f.write('       | {0[7]:7.1f} {0[8]:7.1f} {0[9]:7.1f} {0[10]:7.2f} {0[11]:7.3f} {0[12]:7.3f} {0[13]:7.3f} |\n'.format(T))
                f.write('  +/-  | {0[7]:7.1f} {0[8]:7.1f} {0[9]:7.1f} {0[10]:7.2f} {0[11]:7.3f} {0[12]:7.3f} {0[13]:7.3f} |\n'.format(sT))
                f.write('-------|---------------------------------------------------------|\n')
            f.write('\n')

            # Write residuals
            f.write('Residuals :\n')
            f.write('-----------\n')
            f.write('\n')

            # 1st case : 14 parameters
            if (with_vel):
                f.write('             |         Raw residuals (mm, mm/y)          |           Normalized residuals            |\n')
                f.write('-------------|-------------------------------------------|-------------------------------------------|\n')
                f.write('code pt soln |    E      N      H     VE     VN     VH   |    E      N      H     VE     VN     VH   |\n')
                f.write('-------------|-------------------------------------------|-------------------------------------------|\n')
                for i in range(nsta):
                    f.write('{0} {1} {2} |'.format(code[i], pt[i], soln[i]))
                    f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} {0[3]:6.1f} {0[4]:6.1f} {0[5]:6.1f} |'.format(vl[i]))
                    f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} {0[3]:6.1f} {0[4]:6.1f} {0[5]:6.1f} |\n'.format(vn[i]))
                f.write('-------------|-------------------------------------------|-------------------------------------------|\n')

            # 2nd case : 7 parameters
            else:
                f.write('             |  Raw residuals (mm)  | Normalized residuals |\n')
                f.write('-------------|----------------------|----------------------|\n')
                f.write('code pt soln |    E      N      H   |    E      N      H   |\n')
                f.write('-------------|----------------------|----------------------|\n')
                for i in range(nsta):
                    f.write('{0} {1} {2} |'.format(code[i], pt[i], soln[i]))
                    f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} |'.format(vl[i]))
                    f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} |\n'.format(vn[i]))
                f.write('-------------|----------------------|----------------------|\n')
            f.write('\n')

            f.close()

        return (T, Q, code, pt, soln, vl, vn, w)

#-------------------------------------------------------------------------------
# Routine : trim_metadata
# Purpose : Remove metadata that are not relevant for specified period
# Author  : P. Rebischung
# Created : 18-May-2012
#
# Changes :
#
# Input   : - start : Start of period of interest in SINEX format
#           - end   : End of period of interest in SINEX format
# Output  :
#-------------------------------------------------------------------------------
    def trim_metadata(self, start, end):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Loop over stations
        for s in snx.sta:

            # Loop over receivers
            i = 0
            while (i < len(s.receiver)):
                r = s.receiver[i]

                # If current receiver is relevant for specified period, keep it.
                if (((earlier(r.start, end)) or (r.start == '00:000:00000')) and
                    ((earlier(start, r.end)) or (r.end   == '00:000:00000'))):
                    i = i+1

                # Else, remove it.
                else:
                    s.receiver.pop(i)

            # Loop over antennas
            i = 0
            while (i < len(s.antenna)):
                r = s.antenna[i]

                # If current antenna is relevant for specified period, keep it.
                if (((earlier(r.start, end)) or (r.start == '00:000:00000')) and
                    ((earlier(start, r.end)) or (r.end   == '00:000:00000'))):
                    i = i+1

                # Else, remove it.
                else:
                    s.antenna.pop(i)

            # Loop over eccentricities
            i = 0
            while (i < len(s.eccentricity)):
                r = s.eccentricity[i]

                # If current eccentricity is relevant for specified period, keep it.
                if (((earlier(r.start, end)) or (r.start == '00:000:00000')) and
                    ((earlier(start, r.end)) or (r.end   == '00:000:00000'))):
                    i = i+1

                # Else, remove it.
                else:
                    s.eccentricity.pop(i)

#-------------------------------------------------------------------------------
# Routine : update_stalist
# Purpose : Given an AC solution, update list of stations in the IGS combined solution
# Author  : P. Rebischung
# Created : 16-May-2012
#
# Changes :
#
# Input   : - sta : Station list
#           - ac  : AC combination options
# Output  :
#-------------------------------------------------------------------------------
    def update_stalist(self, sta, ac):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Loop over stations in AC solution
        for s in snx.sta:

            # If current station is not yet in station list,
            if (not(s.code+s.pt in [r.code+r.pt for r in sta])):

                r = record()
                r.code = s.code
                r.pt = s.pt
                r.domes = s.domes
                r.tech = s.tech
                r.lon = s.lon
                r.lat = s.lat
                r.h = s.h
                r.description = s.description
                r.soln = []
                r.X = snx.get_xyz([s.code], [s.pt])[0]
                r.input = []

                i = record()
                i.ac = ac.name
                i.receiver = s.receiver
                i.antenna = s.antenna
                i.eccentricity = s.eccentricity
                r.input.append(i)

                sta.append(r)

            # Else (current station already in station list),
            else:
                ind = [r.code+r.pt for r in sta].index(s.code+s.pt)

                i = record()
                i.ac = ac.name
                i.receiver = s.receiver
                i.antenna = s.antenna
                i.eccentricity = s.eccentricity
                sta[ind].input.append(i)

#-------------------------------------------------------------------------------
# Routine : get_core
# Purpose : Get list of core stations in a solution
# Author  : P. Rebischung
# Created : 15-Nov-2013
#
# Changes : PR, 13-Feb-2014 : Add possibility to reject stations with abnormally large formal errors
#
# Input   : - file   : File containing list of core stations (IGb08_core.txt)
#           - ref    : Datum (sinex object)
#           - thr    : Threshold for rejecting stations whose formal errors exceed
#                      (thr * median of formal errors) either in East, North or Up.
#                      Default is None (no rejection).
#           - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def get_core(self, file, ref, thr=None, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Read list of core stations
        core = []
        f = open(file)
        line = f.readline()
        while (line):
            core.append(line.strip().split())
            line = f.readline()
        f.close()

        # Open log file
        if (log):
            if (append):
                f = open(log, 'a')
            else:
                f = open(log, 'w')

            # Write header
            f.write('================================================================================\n')
            f.write('sinex::get_core : Get list of core stations in a solution\n')
            f.write('================================================================================\n')
            f.write('\n')

        # Build list of core stations available in snx : loop over core clusters
        code = []
        pt   = []
        soln = []
        for i in range(len(core)):
            j = 0
            k = 0
            b = False

            # Look for any station of current cluster in snx and ref
            while ((not(b)) and (j < len(core[i]))):
                if (core[i][j] in [p.code for p in snx.param]):
                    isnx = [p.type+p.code for p in snx.param].index('STAX  '+core[i][j])
                    if ('STAX  '+core[i][j]+snx.param[isnx].pt+snx.param[isnx].soln in [p.type+p.code+p.pt+p.soln for p in ref.param]):
                        b = True
                    else:
                        j = j+1
                else:
                    j = j+1

            # If one station of the cluster was found in snx and ref,
            if (b):
                code.append(core[i][j])
                pt.append(snx.param[isnx].pt)
                soln.append(snx.param[isnx].soln)
                if (log):
                    f.write('{0} : core station found\n'.format(core[i][j]))

            # Else,
            elif (log):
                f.write('{0} : not available\n'.format(core[i][0]))

        # If requested, reject stations with abnormally large formal errors
        if (thr):

            # Compute 3D formal errors
            sig = numpy.zeros(len(code))
            for i in range(len(code)):
                ind = [p.type+p.code+p.pt+p.soln for p in snx.param].index('STAX  '+code[i]+pt[i]+soln[i])
                sig[i] = sqrt(snx.sig[ind+0]**2 + snx.sig[ind+1]**2 + snx.sig[ind+2]**2)

            # Get indices of outliers
            ind = numpy.nonzero(sig > thr*numpy.median(sig))[0]

            # If any outlier,
            if (len(ind) > 0):

                # Write outliers into log file
                if (log):
                    f.write('\n')
                    for i in range(len(ind)):
                        f.write('{0} : Rejected because of abnormally large formal error\n'.format(code[ind[i]]))

                # Reject outliers
                indk = numpy.setdiff1d(list(range(len(code))), ind)
                code = [code[i] for i in indk]
                pt   = [pt[i] for i in indk]
                soln = [soln[i] for i in indk]

        # Close log file
        if (log):
            f.write('\n')
            f.close()

        return (code, pt, soln)

#-------------------------------------------------------------------------------
# Routine : extract_core
# Purpose : Extract core stations from a solution
# Author  : P. Rebischung
# Created : 15-Nov-2013
#
# Changes : PR, 13-Feb-2014 : Add possibility to reject stations with abnormally large formal errors
#
# Input   : - file   : File containing list of core stations (IGb08_core.txt)
#           - ref    : Datum (sinex object)
#           - thr    : Threshold for rejecting stations whose formal errors exceed
#                      (thr * median of formal errors) either in East, North or Up.
#                      Default is None (no rejection).
#           - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def extract_core(self, file, ref, thr=None, log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Get list of core stations
        (code, pt, soln) = snx.get_core(file, ref, thr, log, append)

        # Remove non-core stations from snx
        snx.keep_sta(code, pt, soln, log, append)

#-------------------------------------------------------------------------------
# Routine : calib_lod
# Purpose : Calibrate LOD estimates wrt reference series (Bulletin A)
# Author  : P. Rebischung
# Created : 17-Dec-2011
#
# Changes : PR, 25-Sep-2014 : Apply correction at the normal equation level
#
# Input   : - rec    : erp object or file containing time series of LOD estimates
#           - ref    : Reference erp object or file (Bulletin A)
#           - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
#     def calib_lod(self, rec, ref, log=None, append=False):
# 
#         # snx = self for class methof # Added JMN 05/01/2018
#         snx = self
# 
# 
#         # If necessary, open log file and write header
#         if (log):
#             if (append):
#                 f = open(log, 'a')
#             else:
#                 f = open(log, 'w')
# 
#             # Write header
#             f.write('================================================================================\n')
#             f.write('sinex::calib_lod : Calibrate LOD estimates wrt reference series (Bulletin A)\n')
#             f.write('================================================================================\n')
#             f.write('\n')
# 
#         # Read ERP time series
#         if (rec.__class__.__name__ != 'erp'):
#             rec = erp.read_igs(rec)
# 
#         # Read reference ERP time series
#         if (ref.__class__.__name__ != 'erp'):
#             ref = erp.read_bulla(ref)
# 
#         # Loop over LOD parameters in sinex object
#         for i in range(snx.npar):
#             p = snx.param[i]
#             if (p.type == 'LOD   '):
# 
#                 # Initializations
#                 b = 0
#                 n = 0
#                 mjdref = (date.from_snxepoch(p.tref)).mjd
# 
#                 # Loop over the 10 previous days
#                 for mjd in numpy.arange(mjdref-10, mjdref):
#                     if ((mjd in rec.mjd) and (mjd in ref.mjd)):
#                         irec = (numpy.nonzero(rec.mjd == mjd))[0][0]
#                         iref = (numpy.nonzero(ref.mjd == mjd))[0][0]
#                         b = b + ref.lod[iref] - rec.lod[irec]
#                         n = n+1
# 
#                 # If at least one previous day was available,
#                 if (n > 0):
# 
#                     # Modify normal equation
#                     b = b / n
#                     snx.k = snx.k + b * snx.N[:,i]
# 
#                     # Write message in log file
#                     if (log):
#                         f.write('LOD    {0.soln} {0.tref} corrected by {1:21.14e} ms (mean over {2:2d} days)\n'.format(p, b, n))
# 
#                 # Else, write message in log file
#                 elif (log):
#                     f.write('LOD    {0.soln} {0.tref} NOT corrected !\n'.format(p))
# 
#         # Close log file if necessary
#         if (log):
#             f.write('\n')
#             f.close()

#-------------------------------------------------------------------------------
# Routine : apriori_sf
# Purpose : Compute a priori scale factor of a solution
# Author  : P. Rebischung
# Created : 16-May-2012
#
# Changes :
#
# Input   : - stdref : Reference 3D formal error (mm)
# Output  : - sf     : A priori scale factor
#-------------------------------------------------------------------------------
    def apriori_sf(self, stdref):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Get list of points in solution
        ind  = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'STAX  ')[0]
        code = [snx.param[i].code for i in ind]
        pt   = [snx.param[i].pt for i in ind]
        soln = [snx.param[i].soln for i in ind]

        # Get 3D formal errors
        senh = snx.get_sigenh(code, pt, soln)
        s3d  = numpy.sqrt(numpy.sum(senh**2, axis=1))

        # Compute a priori scale factor, given reference 3D formal error
        sf = stdref / numpy.median(s3d) / 1000

        return sf

#-------------------------------------------------------------------------------
# Routine : add_metadata
# Purpose : Add metadata blocks into a sinex object
# Author  : P. Rebischung
# Created : 09-Mar-2012
#
# Changes :
#
# Input   : - metasnx  : sinex object with information from sitelogs
#           - ref      : File with information for FILE/REFERENCE block. Default is None.
#           - comments : File with information for FILE/COMMENTS block. Default is None.
#           - acks     : File with information for INPUT/ACKNOWLEDGEMENTS block. Default is None.
#           - t        : date object (needed for FILE/REFERENCE block). Default is None.
# Output  :
#-------------------------------------------------------------------------------
    def add_metadata(self, metasnx, ref=None, comments=None, acks=None, t=None):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Fill FILE/REFERENCE block if requested
        if (ref):
            snx.ref = read_yaml(ref, sed=True, t=t)

            # Replace agency if specified in the reference record
            if (snx.ref.agency):
                snx.agency = snx.ref.agency

        # Fill FILE/COMMENTS block if requested
        if (comments):
            snx.comment = read_yaml(comments)

        # Fill INPUT/ACKNOWLEDGEMENTS if requested
        if (acks):
            snx.acks = read_yaml(acks)

        # Loop over stations with available metadata
        for s in snx.sta:
            if (s.code in [metas.code for metas in metasnx.sta]):

                # Index of current station in metasnx.sta
                im = [metas.code for metas in metasnx.sta].index(s.code)
                metas = metasnx.sta[im]

                # Copy station metadata from metasnx to snx
                s.receiver     = metas.receiver
                s.antenna      = metas.antenna
                s.eccentricity = metas.eccentricity

        ## Loop over antennas to fill snx.antenna
        #snx.antenna = []
        #for s in snx.sta:
            #for a in s.antenna:

                ## If current antenna is not yet in snx.antenna,
                #if (not(a.type in [ant.type for ant in snx.antenna])):

                    ## And if current antenna is in metasnx.antenna,
                    #if (a.type in [ant.type for ant in metasnx.antenna]):

                        ## Add current antenna to snx.antenna
                        #ind = [ant.type for ant in metasnx.antenna].index(a.type)
                        #snx.antenna.append(metasnx.antenna[ind])

        # Copy metasnx.antenna to snx.antenna
        snx.antenna = metasnx.antenna

#-------------------------------------------------------------------------------
# Routine : get_residuals
# Purpose : Compare solution to a reference solution and get necessary information for res file
# Author  : P. Rebischung
# Created : 20-May-2012
#
# Changes :
#
# Input   : - ref  : Reference solution (sinex object)
#           - erp  : True if ERPs should be compared. Default is False.
#           - gc   : True if origins should be compared. Default is False.
# Output  : - code : 4-char codes of common points
#           - pt   : PT codes of common points
#           - soln : Solution numbers of common points
#           - v    : E,N,H station position residuals (mm)
#           - s    : E,N,H station position formal errors (mm)
#           - verp : ERP residuals
#           - serp : ERP formal errors
#           - vgc  : Geocenter residuals (mm)
#           - sgc  : Geocenter formal errors (mm)
#-------------------------------------------------------------------------------
    def get_residuals(self, ref, erp=False, gc=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Helmert comparison
        (T, Q, code, pt, soln, v, vn, w) = snx.helmert_wrt(ref, params='RT', weighting='full')

        # Get E,N,H station position formal errors
        s = 1000 * snx.get_sigenh(code, pt, soln)

        # If ERPs should be compared,
        verp = None
        serp = None
        if (erp):

            # Initializations
            nerp = len(numpy.nonzero(numpy.array([p.type for p in ref.param]) == 'XPO   ')[0])
            verp = numpy.inf * numpy.ones((6, nerp))
            serp = numpy.inf * numpy.ones((6, nerp))
            ixpo  = -1
            ixpor = -1
            iypo  = -1
            iypor = -1
            iut   = -1
            ilod  = -1

            # Loop over ERPs
            for i in range(len(snx.param)):
                p = snx.param[i]
                if (p.type in ['XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 'UT    ', 'LOD   ']):

                    # Get index of corresponding parameter in reference solution
                    ind = [q.type+q.tref for q in ref.param].index(p.type+p.tref)

                    # Store corresponding residual and formal error
                    if (p.type == 'XPO   '):
                        ixpo = ixpo + 1
                        verp[0,ixpo] = snx.x[i] - ref.x[ind] - T[5]
                        serp[0,ixpo] = snx.sig[i]
                    elif (p.type == 'XPOR  '):
                        ixpor = ixpor + 1
                        verp[1,ixpor] = snx.x[i] - ref.x[ind]
                        serp[1,ixpor] = snx.sig[i]
                    elif (p.type == 'YPO   '):
                        iypo = iypo + 1
                        verp[2,iypo] = snx.x[i] - ref.x[ind] - T[4]
                        serp[2,iypo] = snx.sig[i]
                    elif (p.type == 'YPOR  '):
                        iypor = iypor + 1
                        verp[3,iypor] = snx.x[i] - ref.x[ind]
                        serp[3,iypor] = snx.sig[i]
                    elif (p.type == 'UT    '):
                        iut = iut + 1
                        verp[4,iut] = snx.x[i] - ref.x[ind] + T[6]/(15*dera_dt)
                        serp[4,iut] = snx.sig[i]
                    elif (p.type == 'LOD   '):
                        ilod = ilod + 1
                        verp[5,ilod] = snx.x[i] - ref.x[ind]
                        serp[5,ilod] = snx.sig[i]

        # If geocenters should be compared,
        vgc = None
        sgc = None
        if (gc):

            # If solution contains an XGC parameter
            if ('XGC   ' in [p.type for p in snx.param]):

                # Get indices of XGC parameters in both solutions
                ind  = [p.type for p in snx.param].index('XGC   ')
                indr = [p.type for p in ref.param].index('XGC   ')

                # Geocenter residuals
                vgc = 1000 * (snx.x[ind:ind+3] - ref.x[indr:indr+3]) - T[0:3]

                # Geocenter formal errors
                sgc = 1000 * snx.sig[ind:ind+3]

        return (code, pt, soln, v, s, verp, serp, vgc, sgc)

#-------------------------------------------------------------------------------
# Routine : check_sat_pco
# Purpose : Check satellite antenna phase center offsets
# Author  : P. Rebischung
# Created : 18-Nov-2014
#
# Changes :
#
# Input   : - ref  : Reference sinex object containing nominal satellite PCOs
# Output  :
#-------------------------------------------------------------------------------
    def check_sat_pco(self, ref):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Loop over a priori SATA parameters
        for i in range(len(snx.prior)):
            p = snx.prior[i]
            if (p.type[0:4] == 'SATA'):

                # Index of corresponding parameter in reference solution
                ind = [pref.type+pref.code for pref in ref.param].index(p.type+p.code)

                # If a priori value differs from reference value,
                if (snx.x0[i] != ref.x[ind]):
                    print('{0} {1} : {2:21.14e} instead of {3:21.14e}'.format(p.type, p.code, snx.x0[i], ref.x[ind]))

#-------------------------------------------------------------------------------
# Routine : refsys_effect
# Purpose : Compute generalized reference system effect of a solution
# Author  : P. Rebischung
# Created : 17-Nov-2011
#
# Changes :
#
# Input   : - params : List of keywords indicating on what "parameters" the reference
#                      system effect should be computed. It can contain the following
#                      keywords:
#                      - 'T'    (How well is defined the origin of the solution?)
#                      - 'S'    (How well is defined the scale of the solution?)
#                      - 'R'    (How well is defined the orientation of the solution?)
#                      - 'dT'   (How well is defined the origin rate of the solution?)
#                      - 'dS'   (How well is defined the scale rate of the solution?)
#                      - 'dR'   (How well is defined the orientation rate of the solution?)
#                      - 'XPO'  (How well is defined the mean of the XPO parameters?)
#                      - 'XPOR' (How well is defined the mean of the XPOR parameters?)
#                      - 'YPO'  (How well is defined the mean of the YPO parameters?)
#                      - 'YPOR' (How well is defined the mean of the YPOR parameters?)
#                      - 'UT'   (How well is defined the mean of the UT parameters?)
#                      - 'LOD'  (How well is defined the mean of the LOD parameters?)
#                      - 'XPCO' (How well is defined the mean of the SATA_X parameters?)
#                      - 'YPCO' (How well is defined the mean of the SATA_Y parameters?)
#                      - 'ZPCO' (How well is defined the mean of the SATA_Z parameters?)
#                      - 'GC'   (How well are defined the XGC/YGC/ZGC parameters?)
#                      - 'mT'   (How well are defined the means of the TX/TY/TZ parameters - in case of a combined solution?)
#                      - 'mS'   (How well is defined the mean of the SC parameters - in case of a combined solution?)
#                      - 'mR'   (How well are defined the means of the RX/RY/RZ parameters - in case of a combined solution?)
#                      - 'all'  (Compute reference system effect on all possible "parameters". This is the default.)
#           - log    : Log file. Default is None.
#           - append : If true, text will be appended to log file. Default if False.
# Output  :
#-------------------------------------------------------------------------------
    def refsys_effect(self, params=['all'], log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initializations
        par_names = ['TX', 'TY', 'TZ', 'SC', 'RX', 'RY', 'RZ', 'dTX', 'dTY', 'dTZ', 'dSC', 'dRX', 'dRY', 'dRZ',
                     'mXPO', 'mXPOR', 'mYPO', 'mYPOR', 'mUT', 'mLOD', 'mXPCO', 'mYPCO', 'mZPCO',
                     'XGC', 'YGC', 'ZGC', 'mTX', 'mTY', 'mTZ', 'mSC', 'mRX', 'mRY', 'mRZ']
        par_units = ['mm', 'mm', 'mm', 'ppb', 'mas', 'mas', 'mas', 'mm/y', 'mm/y', 'mm/y', 'ppb/y', 'mas/y', 'mas/y', 'mas/y',
                     'mas', 'mas/d', 'mas', 'mas/d', 'ms', 'ms/d', 'mm', 'mm', 'mm',
                     'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'ppb', 'mas', 'mas', 'mas']
        x = None
        if (type(snx.x).__name__ != 'NoneType'):
            x = snx.x
        else:
            x = snx.x0
        A = numpy.zeros((len(snx.param), len(par_names)))

        # Fill design matrix
        for i in range(len(snx.param)):

            if (snx.param[i].type == 'STAX  '):
                A[i+0, 0] =  1e-3
                A[i+1, 1] =  1e-3
                A[i+2, 2] =  1e-3
                A[i+0, 3] =  1e-9 * x[i+0]
                A[i+1, 3] =  1e-9 * x[i+1]
                A[i+2, 3] =  1e-9 * x[i+2]
                A[i+1, 4] = -mas2rad * x[i+2]
                A[i+2, 4] =  mas2rad * x[i+1]
                A[i+0, 5] =  mas2rad * x[i+2]
                A[i+2, 5] = -mas2rad * x[i+0]
                A[i+0, 6] = -mas2rad * x[i+1]
                A[i+1, 6] =  mas2rad * x[i+0]

            elif (snx.param[i].type == 'VELX  '):
                A[i+0, 7]  =  1e-3
                A[i+1, 8]  =  1e-3
                A[i+2, 9]  =  1e-3
                A[i+0, 10] =  1e-9 * x[i-3]
                A[i+1, 10] =  1e-9 * x[i-2]
                A[i+2, 10] =  1e-9 * x[i-1]
                A[i+1, 11] = -mas2rad * x[i-1]
                A[i+2, 11] =  mas2rad * x[i-2]
                A[i+0, 12] =  mas2rad * x[i-1]
                A[i+2, 12] = -mas2rad * x[i-3]
                A[i+0, 13] = -mas2rad * x[i-2]
                A[i+1, 13] =  mas2rad * x[i-3]

            elif (snx.param[i].type == 'XPO   '):
                A[i+0, 14] = 1
            elif (snx.param[i].type == 'XPOR  '):
                A[i+0, 15] = 1
            elif (snx.param[i].type == 'YPO   '):
                A[i+0, 16] = 1
            elif (snx.param[i].type == 'YPOR  '):
                A[i+0, 17] = 1
            elif (snx.param[i].type == 'UT    '):
                A[i+0, 18] = 1
            elif (snx.param[i].type == 'LOD   '):
                A[i+0, 19] = 1

            elif (snx.param[i].type == 'SATA_X'):
                A[i+0, 20] = 1e-3
            elif (snx.param[i].type == 'SATA_Y'):
                A[i+0, 21] = 1e-3
            elif (snx.param[i].type == 'SATA_Z'):
                A[i+0, 22] = 1e-3

            elif (snx.param[i].type == 'XGC   '):
                A[i+0, 23] = 1e-3
            elif (snx.param[i].type == 'YGC   '):
                A[i+0, 24] = 1e-3
            elif (snx.param[i].type == 'ZGC   '):
                A[i+0, 25] = 1e-3

            elif (snx.param[i].type == 'TX    '):
                A[i+0, 26] = 1e-3
            elif (snx.param[i].type == 'TY    '):
                A[i+0, 27] = 1e-3
            elif (snx.param[i].type == 'TZ    '):
                A[i+0, 28] = 1e-3
            elif (snx.param[i].type == 'SC    '):
                A[i+0, 29] = 1
            elif (snx.param[i].type == 'RX    '):
                A[i+0, 30] = 1
            elif (snx.param[i].type == 'RY    '):
                A[i+0, 31] = 1
            elif (snx.param[i].type == 'RZ    '):
                A[i+0, 32] = 1

        # Get indices of required parameters
        ind = []

        if ('T' in params):
            ind.extend(list(range(0, 3)))
        if ('S' in params):
            ind.append(3)
        if ('R' in params):
            ind.extend(list(range(4, 7)))

        if ('dT' in params):
            ind.extend(list(range(7, 10)))
        if ('dS' in params):
            ind.append(10)
        if ('dR' in params):
            ind.extend(list(range(11, 14)))

        if ('XPO' in params):
            ind.append(14)
        if ('XPOR' in params):
            ind.append(15)
        if ('YPO' in params):
            ind.append(16)
        if ('YPOR' in params):
            ind.append(17)
        if ('UT' in params):
            ind.append(18)
        if ('LOD' in params):
            ind.append(19)

        if ('XPCO' in params):
            ind.append(20)
        if ('YPCO' in params):
            ind.append(21)
        if ('ZPCO' in params):
            ind.append(22)

        if ('GC' in params):
            ind.extend(list(range(23, 26)))

        if ('mT' in params):
            ind.extend(list(range(26, 29)))
        if ('mS' in params):
            ind.append(29)
        if ('mR' in params):
            ind.extend(list(range(30, 33)))

        # Restrict design matrix to required parameters
        A = A[:,ind]

        # Solution's normal matrix
        N = None
        if (type(snx.N).__name__ != 'NoneType'):
            N = snx.N
        else:
            N = syminv(snx.Q)

        # Covariance matrix of reference frame parameters
        Q = dot(dot(A.T, N), A)
        Q = syminv(Q)

        ## Covariance matrix of reference frame parameters
        #B = dot(syminv(dot(A.T, A)), A.T)
        #Q = dot(B, dot(snx.Q, B.T))

        # Correlation matrix
        C = Q2C(Q)

        # Put correlations in a vector
        name1 = []
        name2 = []
        corr  = []
        for i in range(len(Q)):
            for j in range(i):
                name1.append(par_names[ind[i]])
                name2.append(par_names[ind[j]])
                corr.append(C[i,j])
        corr = numpy.array(corr)

        # Sort correlations by absolute values
        ic = numpy.argsort(numpy.abs(corr))

        # Print results - 1st case : No log file
        if (log is None):

            print('')
            print('Standard deviations of meta-parameters:')
            print('---------------------------------------')
            print('')
            for i in range(len(ind)):
                if (par_units[ind[i]] == 'mm'):
                    print('{0:>5s} : {1:15.8f}    mm ~ {1:15.8f}    mm'.format(par_names[ind[i]], sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'ppb'):
                    print('{0:>5s} : {1:15.8f}   ppb ~ {2:15.8f}    mm'.format(par_names[ind[i]], sqrt(Q[i,i]), 1e-6*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'mas'):
                    print('{0:>5s} : {1:15.8f}   mas ~ {2:15.8f}    mm'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*mas2rad*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'mm/y'):
                    print('{0:>5s} : {1:15.8f}  mm/y ~ {1:15.8f}  mm/y'.format(par_names[ind[i]], sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'ppb/y'):
                    print('{0:>5s} : {1:15.8f} ppb/y ~ {2:15.8f}  mm/y'.format(par_names[ind[i]], sqrt(Q[i,i]), 1e-6*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'mas/y'):
                    print('{0:>5s} : {1:15.8f} mas/y ~ {2:15.8f}  mm/y'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*mas2rad*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'mas/d'):
                    print('{0:>5s} : {1:15.8f} mas/d ~ {2:15.8f}  mm/d'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*mas2rad*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'ms'):
                    print('{0:>5s} : {1:15.8f}    ms ~ {2:15.8f}    mm'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*ms2rad*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'ms/d'):
                    print('{0:>5s} : {1:15.8f}  ms/d ~ {2:15.8f}  mm/d'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*ms2rad*ae*sqrt(Q[i,i])))
            print('')

            print('Correlations:')
            print('-------------')
            print('')
            for i in range(len(corr)-1, -1, -1):
                print('{0:>5s} - {1:>5s} : {2:13.10f}'.format(name1[ic[i]], name2[ic[i]], corr[ic[i]]))
            print('')

        # Print results - 2nd case : Log file
        else:
            if (append):
                f = open(log, 'a')
            else:
                f = open(log, 'w')

            f.write('================================================================================\n')
            f.write('sinex::refsys_effect : Compute generalized reference system effect of a solution\n')
            f.write('================================================================================\n')
            f.write('\n')

            f.write('Standard deviations of meta-parameters:\n')

            f.write('---------------------------------------\n')
            f.write('\n')
            for i in range(len(ind)):
                if (par_units[ind[i]] == 'mm'):
                    f.write('{0:>5s} : {1:15.8f}    mm ~ {1:15.8f}    mm\n'.format(par_names[ind[i]], sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'ppb'):
                    f.write('{0:>5s} : {1:15.8f}   ppb ~ {2:15.8f}    mm\n'.format(par_names[ind[i]], sqrt(Q[i,i]), 1e-6*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'mas'):
                    f.write('{0:>5s} : {1:15.8f}   mas ~ {2:15.8f}    mm\n'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*mas2rad*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'mm/y'):
                    f.write('{0:>5s} : {1:15.8f}  mm/y ~ {1:15.8f}  mm/y\n'.format(par_names[ind[i]], sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'ppb/y'):
                    f.write('{0:>5s} : {1:15.8f} ppb/y ~ {2:15.8f}  mm/y\n'.format(par_names[ind[i]], sqrt(Q[i,i]), 1e-6*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'mas/y'):
                    f.write('{0:>5s} : {1:15.8f} mas/y ~ {2:15.8f}  mm/y\n'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*mas2rad*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'mas/d'):
                    f.write('{0:>5s} : {1:15.8f} mas/d ~ {2:15.8f}  mm/d\n'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*mas2rad*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'ms'):
                    f.write('{0:>5s} : {1:15.8f}    ms ~ {2:15.8f}    mm\n'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*ms2rad*ae*sqrt(Q[i,i])))
                elif (par_units[ind[i]] == 'ms/d'):
                    f.write('{0:>5s} : {1:15.8f}  ms/d ~ {2:15.8f}  mm/d\n'.format(par_names[ind[i]], sqrt(Q[i,i]), 1000*ms2rad*ae*sqrt(Q[i,i])))
            f.write('\n')

            f.write('Correlations:\n')
            f.write('-------------\n')
            f.write('\n')
            for i in range(len(corr)-1, -1, -1):
                f.write('{0:>5s} - {1:>5s} : {2:13.10f}\n'.format(name1[ic[i]], name2[ic[i]], corr[ic[i]]))
            f.write('\n')

#-------------------------------------------------------------------------------
# Routine : correct_opoleload
# Purpose : Correct ocean pole tide loading displacements in a normal equation
# Author  : P. Rebischung
# Created : 10-Feb-2015
#
# Changes :
#
# Input   : - opole   : Object containing interpolating splines for the 6 OPTL coefficients
#           - mjd     : MJD
#           - meanpm  : Keyword indicating which mean pole model should be used
#                       ('IERS2003' or 'IERS2010'). Default is 'IERS2010'.
#           - reverse : True if ocean pole tide displacements should be added back instead
#                       of removed. Default is False.
# Output  :
#-------------------------------------------------------------------------------
    def correct_opoleload(self, opole, mjd, meanpm='IERS2010', reverse=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Get indices of STAX parameters
        indx = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'STAX  ')[0]

        # Get station coordinates and XYZ->ENH rotation matrices
        X = numpy.zeros((len(indx), 3))
        for i in range(len(indx)):
            X[i,:] = snx.x[indx[i]:indx[i]+3]
        (phi, lam, h) = cart2geo(X)
        lat = 180/pi*phi
        lon = 180/pi*lam
        R = xyz2enh(X)

        # Get pole coordinates
        ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'XPO   ')[0][0]
        xpo = snx.x[ind]
        ind = numpy.nonzero(numpy.array([p.type for p in snx.param]) == 'YPO   ')[0][0]
        ypo = snx.x[ind]

        # Compute ENH ocean pole tide loading displacements
        denh = compute_opoleload(opole, lon, lat, mjd, xpo, ypo, meanpm)

        # Compute vector of XYZ displacements
        dx = numpy.zeros(snx.npar)
        for i in range(len(indx)):
            dx[indx[i]:indx[i]+3] = dot(R[i].T, denh[i,:])

        # If necessary, change sign of displacements
        if (reverse):
            dx = -dx

        # Modify 2nd member of normal equation
        snx.k = snx.k - dot(snx.N, dx)

        ## Modify estimated station positions
        #snx.x = snx.x - dx

#-------------------------------------------------------------------------------
# Routine : get_nonpubsolns
# Purpose : Get list of solns to remove from official IGS cumulative solution
# Author  : P. Rebischung
# Created : 13-Jun-2012
#
# Changes :
#
# Input   : - sigmax : Maximal velocity formal error (m/y)
#           - fcontr : File with velocity constraints
# Output  : - code   : 4-char codes of solns to remove
#           - pt     : PT codes of solns to remove
#           - soln   : Solution numbers of solns to remove
#           - domes  : DOMES numbers of solns to remove
#-------------------------------------------------------------------------------
    def get_nonpubsolns(self, sigmax, fcontr):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initialization
        code  = []
        pt    = []
        soln  = []
        domes = []

        # Get list of solns with constrained velocities
        f = open(fcontr)
        line = f.readline()
        while (line):
            if (line[28:50] == '0.0000  0.0000  0.0000'):
                code.append(line[64:68])
                pt.append(line[69:71])
                soln.append('  '+line[10:12])
                domes.append(line[0:9])
            line = f.readline()
        f.close()

        # Loop over stations
        for s in snx.sta:
            keep = False

            # Keep current station if it has a DOMES number
            if (s.domes != '     M   '):

                # and at least one velocity with formal error < sigmax
                for p in s.soln:
                    i = [par.type+par.code+par.pt+par.soln for par in snx.param].index('VELX  '+s.code+s.pt+p.soln)
                    sig = sqrt(snx.sig[i+0]**2 + snx.sig[i+1]**2 + snx.sig[i+2]**2)
                    if (sig < sigmax):
                        keep = True

            # If necessary, add all solns of current station to the black list
            if (not(keep)):
                for p in s.soln:
                    code.append(s.code)
                    pt.append(s.pt)
                    soln.append(p.soln)
                    domes.append(s.domes)

        return (code, pt, soln, domes)

#-------------------------------------------------------------------------------
# Routine : add_psd
# Purpose : Add post-seismic deformation models to a solution
# Author  : P. Rebischung
# Created : 10-Nov-2015
#
# Changes :
#
# Input   : psd : sinex object containing post-seismic deformation models
# Output  :
#------------------------------------------------------------------------------
    def add_psd(self, psd):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # List of stations with post-seismic deformation models
        codept = [p.code+p.pt for p in psd.param]

        # Loop over STAX parameters
        for i in range(snx.npar):
            p = snx.param[i]
            if (p.type == 'STAX  '):

                # If current station has post-seismic deformation models,
                if (p.code+p.pt in codept):

                    # Compute ENH post-seismic deformations
                    (denh, senh) = compute_psd(psd, p.code, p.tref)

                    # Compute XYZ post-seismic deformations
                    R = xyz2enh(snx.x[i:i+3])
                    dxyz = dot(R.T, denh)
                    Qenh = numpy.diag(senh**2)
                    Qxyz = dot(R.T, dot(Qenh, R))

                    # Add post-seismic deformations
                    snx.x[i:i+3] = snx.x[i:i+3] + dxyz
                    if (snx.Q != None):
                        snx.Q[i:i+3,i:i+3] = snx.Q[i:i+3,i:i+3] + Qxyz
                        snx.sig[i:i+3] = numpy.sqrt(numpy.diag(snx.Q[i:i+3,i:i+3]))

#-------------------------------------------------------------------------------
# Routine : remove_psd
# Purpose : Remove post-seismic deformation models from a solution
#           (without touching covariance matrix)
# Author  : P. Rebischung
# Created : 16-Jan-2017
#
# Changes :
#
# Input   : psd : sinex object containing post-seismic deformation models
# Output  :
#------------------------------------------------------------------------------
    def remove_psd(self, psd):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # List of stations with post-seismic deformation models
        codept = [p.code+p.pt for p in psd.param]

        # Loop over STAX parameters
        for i in range(snx.npar):
            p = snx.param[i]
            if (p.type == 'STAX  '):

                # If current station has post-seismic deformation models,
                if (p.code+p.pt in codept):

                    # Compute ENH post-seismic deformations
                    (denh, senh) = compute_psd(psd, p.code, p.tref)

                    # Compute XYZ post-seismic deformations
                    R = xyz2enh(snx.x[i:i+3])
                    dxyz = dot(R.T, denh)

                    # Remove post-seismic deformations
                    snx.x[i:i+3] = snx.x[i:i+3] - dxyz

#-------------------------------------------------------------------------------
# Routine : add_per
# Purpose : Add periodic terms to a solution
# Author  : P. Rebischung
# Created : 08-Mar-2016
#
# Changes :
#
# Input   : - file : File containing periodic term coefficients
#           - T    : Period (days)
#           - tref : Reference epoch
# Output  :
#------------------------------------------------------------------------------
    def add_per(self, file, T, tref):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Read input file
        (code, pt, soln, ae, be, an, bn, ah, bh) = numpy.loadtxt(file, dtype='O,O,i,f,f,f,f,f,f', unpack=True)
        pt = [' '+p for p in pt]
        soln = ['{0:4d}'.format(s) for s in soln]
        codeptsoln = [code[i]+pt[i]+soln[i] for i in range(len(code))]

        # Reference MJD
        mjd0 = date.from_snxepoch(tref).mjd

        # Loop over STAX parameters
        for i in range(snx.npar):
            p = snx.param[i]
            if (p.type == 'STAX  '):

                # Get index of corresponding periodic terms
                ind = codeptsoln.index(p.code+p.pt+p.soln)

                # Compute ENH deformations
                dt = date.from_snxepoch(p.tref).mjd - mjd0
                c = cos(2*pi*dt/T)
                s = sin(2*pi*dt/T)
                denh = numpy.zeros(3)
                denh[0] = ae[ind]*c + be[ind]*s
                denh[1] = an[ind]*c + bn[ind]*s
                denh[2] = ah[ind]*c + bh[ind]*s

                # Compute and add XYZ deformations
                R = xyz2enh(snx.x[i:i+3])
                snx.x[i:i+3] = snx.x[i:i+3] + dot(R.T, denh) / 1000

#-------------------------------------------------------------------------------
# Routine : plate_poles
# Purpose : Estimate tectonic plate rotation poles from site velocities
# Author  : P. Rebischung
# Created : 10-Jun-2016
#
# Changes :
#
# Input   : - file      : File containing station <-> plate correspondence
#           - Trate     : True or False depending on whether a translation rate bias
#                         should be estimated. Default is True.
#           - weighting : Keyword to indicate which weighting should be used, i.e.,
#                         'identity', 'diagonal', or 'full' (default).
#           - log       : Log file. Default is None.
#           - append    : If true, text will be appended to log file. Default if False.
# Output  :
#------------------------------------------------------------------------------
    def plate_poles(self, file, Trate=True, weighting='full', log=None, append=False):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Read station <-> plate correspondence file
        code = []
        plate = []
        f = open(file)
        line = f.readline()
        while (line):
            tab = line.strip().split()
            code.append(tab[0])
            plate.append(tab[1])
            line = f.readline()
        f.close()

        # List of plates
        plate_list = list(set(plate))

        # Initializations
        plate_nsta = []
        plate_sta = []
        plate_ivx = []

        # Loop over plates
        i = 0
        while (i < len(plate_list)):

            # Initializations
            sta = []
            ivx = []

            # Loop over VELX parameters
            for j in range(snx.npar):
                p = snx.param[j]
                if (p.type == 'VELX  '):

                    # If current point belongs to current plate,
                    ind = code.index(p.code)
                    if (plate[ind] == plate_list[i]):
                        sta.append(p.code)
                        ivx.append(j)

            # If at least two points on current plate,
            if (len(sta) > 1):

                # Sort stations
                ind = numpy.argsort(sta)
                sta = [sta[k] for k in ind]
                ivx = [ivx[k] for k in ind]

                # Save information
                plate_nsta.append(len(sta))
                plate_sta.append(sta)
                plate_ivx.append(ivx)
                i = i+1

            # Else, remove current plate
            else:
                plate_list.pop(i)

        # Sort plates
        ind = numpy.argsort(plate_list)
        plate_list = [plate_list[i] for i in ind]
        plate_nsta = [plate_nsta[i] for i in ind]
        plate_sta = [plate_sta[i] for i in ind]
        plate_ivx = [plate_ivx[i] for i in ind]

        # Initializations
        nplate = len(plate_list)
        if (Trate):
            npar = 3*nplate + 3
        else:
            npar = 3*nplate
        nsta = sum(plate_nsta)
        iQ = numpy.zeros(3*nsta, dtype='int')
        X = numpy.zeros(3*nsta)
        V = numpy.zeros(3*nsta)
        lon = numpy.zeros(nsta)
        lat = numpy.zeros(nsta)
        A = numpy.zeros((2*nsta, npar))
        R = numpy.zeros((2*nsta, 3*nsta))
        ista = -1

        # Loop over stations
        for i in range(nplate):
            for j in range(plate_nsta[i]):
                ista = ista+1
                iQ[3*ista:3*ista+3] = list(range(plate_ivx[i][j], plate_ivx[i][j]+3))

                # Get station coordinates
                X[3*ista:3*ista+3] = snx.x[plate_ivx[i][j]-3:plate_ivx[i][j]]
                V[3*ista:3*ista+3] = snx.x[plate_ivx[i][j]:plate_ivx[i][j]+3]
                (phi, lam, h) = cart2geo(X[3*ista:3*ista+3])
                lon[ista] = 180/pi*lam
                lat[ista] = 180/pi*phi

                # Update global rotation matrix
                Ri = xyz2enh(X[3*ista:3*ista+3])
                R[2*ista:2*ista+2, 3*ista:3*ista+3] = Ri[0:2,:]

                # VE,VN / pole partial derivatives
                Ai = numpy.zeros((3, 3))
                Ai[1,0] = -X[3*ista+2]
                Ai[2,0] =  X[3*ista+1]
                Ai[0,1] =  X[3*ista+2]
                Ai[2,1] = -X[3*ista+0]
                Ai[0,2] = -X[3*ista+1]
                Ai[1,2] =  X[3*ista+0]
                A[2*ista:2*ista+2, 3*i:3*i+3] = as2rad * dot(Ri[0:2,:], Ai)

                # VE,VN / translation rate partial derivatives
                if (Trate):
                    A[2*ista:2*ista+2, 3*nplate:3*nplate+3] = Ri[0:2,:]

        # Build least-squares system
        y = 1000 * dot(R, V)
        if (weighting == 'full'):
            Q = 1e6 * dot(R, dot(snx.Q[numpy.ix_(iQ, iQ)], R.T))
            P = syminv(Q)
            N = dot(A.T, dot(P, A))
            b = dot(A.T, dot(P, y))
        elif (weighting == 'diagonal'):
            Q = 1e6 * dot(R, dot(snx.Q[numpy.ix_(iQ, iQ)], R.T))
            P = numpy.diag(1./numpy.diag(Q))
            N = dot(A.T, dot(P, A))
            b = dot(A.T, dot(P, y))
        elif (weighting == 'identity'):
            N = dot(A.T, A)
            b = dot(A.T, y)

        # Solve least-squares system
        Qx = syminv(N)
        x = dot(Qx, b)
        vm = dot(A, x)
        v = y - vm
        if (weighting == 'identity'):
            s0 = sqrt(sum(v**2) / (2*nsta-npar))
            Qv = numpy.eye(2*nsta) - dot(A, dot(Qx, A.T))
        else:
            s0 = sqrt(dot(v.T, dot(P, v)) / (2*nsta-npar))
            Qv = Q - dot(A, dot(Qx, A.T))
        Qx = s0**2*Qx
        sx = numpy.sqrt(numpy.diag(Qx))
        Qv = s0**2*Qv
        vn = v / numpy.sqrt(numpy.diag(Qv))

        # Compute East and North WRMS
        i = numpy.arange(0, 2*nsta, 2)
        we = sqrt(numpy.sum(v[i]**2 / numpy.diag(Qv[numpy.ix_(i,i)])) / numpy.sum(1. / numpy.diag(Qv[numpy.ix_(i,i)])))
        i = numpy.arange(1, 2*nsta, 2)
        wn = sqrt(numpy.sum(v[i]**2 / numpy.diag(Qv[numpy.ix_(i,i)])) / numpy.sum(1. / numpy.diag(Qv[numpy.ix_(i,i)])))

        # Write log file
        if (log):
            if (append):
                f = open(log, 'a')
            else:
                f = open(log, 'w')

            # Write header
            f.write('================================================================================\n')
            f.write('sinex::plate_poles : from {0}\n'.format(snx.file))
            f.write('================================================================================\n')
            f.write('\n')

            # Write main options and statistics
            f.write('Number of plates   : {0}\n'.format(nplate))
            f.write('Number of stations : {0}\n'.format(nsta))
            f.write('Translation rate   : {0}\n'.format(Trate))
            f.write('Weighting          : {0}\n'.format(weighting))
            f.write('Variance factor    : {0:.8e}\n'.format(s0))
            f.write('WRMS East          : {0:6.3f} mm/yr\n'.format(we))
            f.write('WRMS North         : {0:6.3f} mm/yr\n'.format(wn))
            f.write('\n')

            # Write translation rate if estimated
            if (Trate):
                f.write('Translation rate:\n')
                f.write('-----------------\n')
                f.write('\n')
                f.write(' TX = {0:7.3f} +/- {1:7.3f} mm/yr\n'.format(x[-3], sx[-3]))
                f.write(' TY = {0:7.3f} +/- {1:7.3f} mm/yr\n'.format(x[-2], sx[-2]))
                f.write(' TZ = {0:7.3f} +/- {1:7.3f} mm/yr\n'.format(x[-1], sx[-1]))
                f.write('\n')

            # Write plate summary
            f.write('Tectonic plates:\n')
            f.write('----------------\n')
            f.write('\n')
            f.write('      | #sta      WRMS E/N   |    OX (mas/yr)        OY (mas/yr)        OZ (mas/yr)   |\n')
            f.write('------|----------------------|--------------------------------------------------------|\n')
            i0 = 0
            for i in range(nplate):
                ind = numpy.arange(i0, i0+2*plate_nsta[i], 2)
                wei = sqrt(numpy.sum(v[ind]**2 / numpy.diag(Qv[numpy.ix_(ind,ind)])) / numpy.sum(1. / numpy.diag(Qv[numpy.ix_(ind,ind)])))
                ind = numpy.arange(i0+1, i0+2*plate_nsta[i], 2)
                wni = sqrt(numpy.sum(v[ind]**2 / numpy.diag(Qv[numpy.ix_(ind,ind)])) / numpy.sum(1. / numpy.diag(Qv[numpy.ix_(ind,ind)])))
                i0 = i0+2*plate_nsta[i]
                f.write(' {0} | {1:4d}   {2:6.3f} {3:6.3f} | {4:6.3f} +/- {5:5.3f}   {6:6.3f} +/- {7:5.3f}   {8:6.3f} +/- {9:5.3f} |\n'.format(plate_list[i], plate_nsta[i], wei, wni, x[3*i+0], sx[3*i+0], x[3*i+1], sx[3*i+1], x[3*i+2], sx[3*i+2]))
            f.write('------|----------------------|--------------------------------------------------------|\n')
            f.write('\n')

            # Write residuals
            f.write('Residuals:\n')
            f.write('----------\n')
            f.write('\n')
            f.write(' code   plat | Raw: E/N, mm/yr | Normalized: E/N |\n')
            f.write('-------------|-----------------|-----------------|\n')
            k = 0
            for i in range(nplate):
                for j in range(plate_nsta[i]):
                    f.write(' {0}   {1} | {2:7.3f} {3:7.3f} | {4:7.3f} {5:7.3f} |\n'.format(plate_sta[i][j], plate_list[i], v[2*k+0], v[2*k+1], vn[2*k+0], vn[2*k+1]))
                    k = k+1
                f.write('-------------|-----------------|-----------------|\n')
            f.write('\n')

            # Close log file
            f.close()

        # Re-arrange outputs
        if (Trate):
            omega = x[:-3].reshape((nplate, 3))
            T = x[-3:]
        else:
            omega = x.reshape((nplate, 3))
            T = numpy.zeros(3)

        code = []
        for i in range(nplate):
            code.extend(plate_sta[i])

        v.resize((nsta, 2))
        vn.resize((nsta, 2))
        vm.resize((nsta, 2))

        # Return
        return (plate_list, omega, T, Qx, code, v, vn, vm)

#-------------------------------------------------------------------------------
# Routine : insert_disc
# Purpose : Insert new "discontinuity" for a given station, by duplicating appropriate soln
# Author  : P. Rebischung
# Created : 23-Jun-2016
#
# Changes :
#
# Input   : - code : 4-char station ID
#           - t    : discontinuity date (SINEX format)
# Output  :
#------------------------------------------------------------------------------
    def insert_disc(self, code, t):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Get station from snx.sta
        ind = [s.code for s in snx.sta].index(code)
        s = snx.sta[ind]

        # Look for which soln to split
        i = 0
        b = False
        while (not(b)):
            if ((earlier(s.soln[i].datastart, t)) and (earlier(t, s.soln[i].dataend))):
                b = True
            else:
                i = i+1

        # Modify end (and mean) date of soln i
        end = s.soln[i].dataend
        s.soln[i].dataend = t
        mjd1 = date.from_snxepoch(s.soln[i].datastart).mjd
        mjd2 = date.from_snxepoch(s.soln[i].dataend).mjd
        s.soln[i].datamean = date.from_mjd((mjd1+mjd2)/2).snxepoch()

        # Insert new soln
        oldsoln = s.soln[i].soln
        newsoln = '{0:4d}'.format(int(s.soln[i].soln)+1)
        r = record()
        r.soln = newsoln
        r.datastart = t
        r.dataend = end
        mjd1 = date.from_snxepoch(r.datastart).mjd
        mjd2 = date.from_snxepoch(r.dataend).mjd
        r.datamean = date.from_mjd((mjd1+mjd2)/2).snxepoch()
        s.soln.insert(i+1, r)

        # Modify solns of subsequent solns
        for j in range(len(s.soln)-1, i+1, -1):
            for p in snx.param:
                if ((p.code == code) and (p.soln == s.soln[j].soln)):
                    p.soln = '{0:4d}'.format(int(s.soln[j].soln)+1)
            s.soln[j].soln = '{0:4d}'.format(int(s.soln[j].soln)+1)

        # Get indices of parameters to be duplicated
        ind = []
        for i in range(snx.npar):
            p = snx.param[i]
            if ((p.code == code) and (p.soln == oldsoln)):
                ind.append(i)

        # Duplicate parameters
        for i in range(len(ind)):
            p = copy.deepcopy(snx.param[ind[i]])
            p.soln = newsoln
            snx.param.append(p)

        # Duplicate parameter values and sigmas
        snx.x = numpy.hstack((snx.x, snx.x[ind]))
        snx.sig = numpy.hstack((snx.sig, snx.sig[ind]))

        # Extend covariance matrix
        if (snx.Q != None):
            snx.Q = numpy.vstack((numpy.hstack((snx.Q, snx.Q[:,ind])), numpy.hstack((snx.Q[ind,:], snx.Q[numpy.ix_(ind,ind)]))))

        # Update number of parameters
        snx.npar = len(snx.param)

#-------------------------------------------------------------------------------
# Routine : apply_offsets
# Purpose : Apply position offsets to specified solns
# Author  : P. Rebischung
# Created : 05-Jul-2016
#
# Changes :
#
# Input   : - code   : 4-char station IDs
#           - pt     : PT codes
#           - soln   : Solns
#           - denh   : E/N/H position offsets (mm)
# Output  :
#------------------------------------------------------------------------------
    def apply_offsets(self, code, pt, soln, denh):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Loop over offsets to apply
        for i in range(len(code)):

            # If current soln not in solution,
            if (not('STAX  '+code[i]+pt[i]+soln[i] in [p.type+p.code+p.pt+p.soln for p in snx.param])):
                print(code[i]+' '+pt[i]+' '+soln[i]+' not found in solution !')

            # Else,
            else:

                # Get index of matching STAX parameter
                ind = [p.type+p.code+p.pt+p.soln for p in snx.param].index('STAX  '+code[i]+pt[i]+soln[i])

                # Apply position offset
                X = snx.x[ind:ind+3]
                R = xyz2enh(X)
                snx.x[ind:ind+3] = snx.x[ind:ind+3] + dot(R.T, denh[i]) / 1000

#-------------------------------------------------------------------------------
# Routine : write_solndomes
# Purpose : Write special soln file with DOMES numbers (needed for IGS combination database)
# Author  : P. Rebischung
# Created : 21-Jul-2015
#
# Changes :
#
# Input   : - solns : Soln table
#           - file  : File to write
# Output  :
#-------------------------------------------------------------------------------
    def write_solndomes(self, solns, file):

        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Open output file
        f = open(file, 'w')

        # Loop over stations
        for sta in snx.sta:

            # If station is in soln table,
            if (sta.code+sta.pt in [s.code+s.pt for s in solns]):

                # Get index of station in soln table
                i = [s.code+s.pt for s in solns].index(sta.code+sta.pt)

                # Write P-solns
                for p in solns[i].P:
                    f.write(' {0.code} {0.pt} {0.domes} {1.soln} P {1.start} {1.end} P - {1.reason}\n'.format(sta, p))

                # Write V-solns
                for v in solns[i].V:
                    f.write(' {0.code} {0.pt} {0.domes} {1.soln} P {1.start} {1.end} V - {1.reason}\n'.format(sta, v))

        # Close output file
        f.close()

#-------------------------------------------------------------------------------
# Routine : loadest_wrt
# Purpose : Estimate surface load coefficients from comparison to a reference solution
# Author  : P. Rebischung
# Created : 31-Jan-2012
#
# Changes :
#
# Input   : - ref       : Solution to which the comparison is made (sinex object)
#           - center    : 'CM', 'CF' or 'CN'
#           - lmax      : Maximal SH degree
#           - lovefile  : File containing load Love numbers
#           - weighting : Keyword to indicate which weighting should be used.
#                         It can take the following values :
#                         - 'identity' to use an identity weight matrix (default)
#                         - 'diagonal' to use a diagonal weight matrix
#                         - 'full' to use a full weight matrix (inv(snx.Q))
#           - log       : Log file. Default is None.
#           - append    : If true, text will be appended to log file. Default if False.
# Output  : - T         : Vector of estimated parameters (rotations + load coefficients)
#           - Q         : Covariance matrix of estimated parameters
#           - code      : 4-char codes of points used for the comparison
#           - pt        : PT codes of points used for the comparison
#           - soln      : solns of points used for the comparison
#           - vl        : E, N, H residuals
#           - vln       : Normalized E, N, H residuals
#           - w         : WRMS of topocentric residuals
#-------------------------------------------------------------------------------
    def loadest_wrt(self, ref, center, lmax, lovefile, weighting='full', log=None, append=False):


        # snx = self for class methof # Added JMN 05/01/2018
        snx = self


        # Initializations
        re   = (ae**2*be)**(1./3.)
        me   = 5.9726e24
        rho  = 3*me / (4*pi*re**3)
        nsta = 0
        ind1 = []
        ind0 = []
        code = []
        pt   = []
        soln = []

        # Read Love numbers
        lovel = []
        loveh = []
        f = open(lovefile)
        line = f.readline()
        while (line):
            if (line[0] != '#'):
                tab = line.split()
                loveh.append(float(tab[1]))
                lovel.append(float(tab[2]))
            line = f.readline()
        f.close()
        lovel = numpy.array(lovel)
        loveh = numpy.array(loveh)

        # Loop over STAX parameters of snx to identify common stations
        for i in range(len(snx.param)):
            p1 = snx.param[i]
            if (p1.type == 'STAX  '):

                # If current station is common to both solutions,
                if ('STAX  '+p1.code+p1.pt+p1.soln in [p0.type+p0.code+p0.pt+p0.soln for p0 in ref.param]):
                    nsta += 1
                    code.append(p1.code)
                    pt.append(p1.pt)
                    soln.append(p1.soln)

                    j = [p0.type+p0.code+p0.pt+p0.soln for p0 in ref.param].index('STAX  '+p1.code+p1.pt+p1.soln)
                    ind1.extend(list(range(i, i+3)))
                    ind0.extend(list(range(j, j+3)))

        # Get reference coordinates of common stations
        nobs = len(ind0)
        X = numpy.zeros((nsta, 3))
        X[:,0] = ref.x[ind0[0:nobs:3]]
        X[:,1] = ref.x[ind0[1:nobs:3]]
        X[:,2] = ref.x[ind0[2:nobs:3]]

        # Compute lat, colat, lon and XYZ->ENH rotation matrices
        (phi, lam, h) = cart2geo(X)
        theta = pi/2 - phi
        R = xyz2enh(X)

        # Compute all spherical harmonics up to degree/order lmax
        C = numpy.zeros((lmax+1, lmax+1, nsta))
        S = numpy.zeros((lmax+1, lmax+1, nsta))
        for l in range(0, lmax+1):
            for m in range(0, l+1):
                Y = special.sph_harm(m, l, lam, theta)
                C[l,m,:] = numpy.real(Y)
                S[l,m,:] = numpy.imag(Y)

        # Initialize design matrix A
        if (center == 'CM'):
            npar = lmax*(lmax+2)+3
        elif (center in ['CF', 'CN']):
            npar = lmax*(lmax+2)+6
        A = numpy.zeros((nobs, npar))

        # XYZ / translations-rotations partial derivatives
        if (center == 'CM'):
            A[1:nobs:3, 0] = -mas2rad * X[:, 2]
            A[2:nobs:3, 0] =  mas2rad * X[:, 1]
            A[0:nobs:3, 1] =  mas2rad * X[:, 2]
            A[2:nobs:3, 1] = -mas2rad * X[:, 0]
            A[0:nobs:3, 2] = -mas2rad * X[:, 1]
            A[1:nobs:3, 2] =  mas2rad * X[:, 0]
        elif (center in ['CF', 'CN']):
            A[0:nobs:3, 0] = 1e-3
            A[1:nobs:3, 1] = 1e-3
            A[2:nobs:3, 2] = 1e-3
            A[1:nobs:3, 3] = -mas2rad * X[:, 2]
            A[2:nobs:3, 3] =  mas2rad * X[:, 1]
            A[0:nobs:3, 4] =  mas2rad * X[:, 2]
            A[2:nobs:3, 4] = -mas2rad * X[:, 0]
            A[0:nobs:3, 5] = -mas2rad * X[:, 1]
            A[1:nobs:3, 5] =  mas2rad * X[:, 0]

        # Loop over degrees
        for l in range(1, lmax+1):

            # Index of load(l,0,C) parameter
            if (center == 'CM'):
                ind = (l-1)*(l+1)+3
            elif (center in ['CF', 'CN']):
                ind = (l-1)*(l+1)+6

            # ENH/load(l,0,C) partial derivatives
            A[0:nobs:3, ind] = 0
            A[1:nobs:3, ind] = -3*lovel[l]/((2*l+1)*rho*numpy.sin(theta))*(l*numpy.cos(theta)*C[l,0,:]-l*sqrt(float(2*l+1)/float(2*l-1))*C[l-1,0,:])
            A[2:nobs:3, ind] =  3*loveh[l]/((2*l+1)*rho)*C[l,0,:]

            # Loop over orders > 0
            for m in range(1, l+1):

                # Index of load(l,m,C) parameter
                if (center == 'CM'):
                    ind = (l-1)*(l+1)+2*m+2
                elif (center in ['CF', 'CN']):
                    ind = (l-1)*(l+1)+2*m+5

                # ENH/load(l,m,C) partial derivatives
                A[0:nobs:3, ind] = -3*m*lovel[l]/((2*l+1)*rho*numpy.sin(theta))*S[l,m,:]
                A[1:nobs:3, ind] = -3*lovel[l]/((2*l+1)*rho*numpy.sin(theta))*(l*numpy.cos(theta)*C[l,m,:]-sqrt(float(2*l+1)/float(2*l-1)*float(l+m)*float(l-m))*C[l-1,m,:])
                A[2:nobs:3, ind] =  3*loveh[l]/((2*l+1)*rho)*C[l,m,:]

                # ENH/load(l,m,S) partial derivatives
                ind = ind+1
                A[0:nobs:3, ind] =  3*m*lovel[l]/((2*l+1)*rho*numpy.sin(theta))*C[l,m,:]
                A[1:nobs:3, ind] = -3*lovel[l]/((2*l+1)*rho*numpy.sin(theta))*(l*numpy.cos(theta)*S[l,m,:]-sqrt(float(2*l+1)/float(2*l-1)*float(l+m)*float(l-m))*S[l-1,m,:])
                A[2:nobs:3, ind] =  3*loveh[l]/((2*l+1)*rho)*S[l,m,:]

        # Rotate ENH/load partial derivatives to obtain XYZ/load partial derivatives
        if (center == 'CM'):
            for i in range(nsta):
                A[3*i+0:3*i+3,3:] = dot(R[i].T, A[3*i+0:3*i+3,3:])
        elif (center in ['CF', 'CN']):
            for i in range(nsta):
                A[3*i+0:3*i+3,6:] = dot(R[i].T, A[3*i+0:3*i+3,6:])

        # 2nd member B
        B = snx.x[ind1] - ref.x[ind0]

        # Compute normal matrix and 2nd member of normal equation
        # Also modify design matrix if we're in CN.
        N  = None
        K  = None
        P  = None
        if (weighting == 'identity'):
            if (center == 'CN'):
                AtP  = A[:,0:6].T
                Q    = syminv(dot(AtP, A[:,0:6]))
                proj = dot(dot(A[:,0:6], Q), AtP)
                A[:,6:] = A[:,6:] - dot(proj, A[:,6:])
            N  = dot(A.T, A)
            K  = dot(A.T, B)
        elif (weighting == 'diagonal'):
            P  = [1 / (snx.param[i].sigma**2) for i in ind1]
            if (center == 'CN'):
                AtP  = A[:,0:6].T * P
                Q    = syminv(dot(AtP, A[:,0:6]))
                proj = dot(dot(A[:,0:6], Q), AtP)
                A[:,6:] = A[:,6:] - dot(proj, A[:,6:])
            N  = A.T * P
            K  = dot(N, B)
            N  = dot(N, A)
        elif (weighting == 'full'):
            P  = syminv(snx.Q[numpy.ix_(ind1,ind1)])
            if (center == 'CN'):
                AtP  = dot(A[:,0:6].T, P)
                Q    = syminv(dot(AtP, A[:,0:6]))
                proj = dot(dot(A[:,0:6], Q), AtP)
                A[:,6:] = A[:,6:] - dot(proj, A[:,6:])
            N  = dot(A.T, P)
            K  = dot(N, B)
            N  = dot(N, A)

        # Solve normal equation
        Q = syminv(N)
        T = dot(Q, K)

        # Compute residuals and variance factor
        vg = B - dot(A, T)
        vPv = None
        sigma0 = None

        if (weighting == 'identity'):
            vPv = dot(vg, vg)
        elif (weighting == 'diagonal'):
            vPv = dot(vg*P, vg)
        elif (weighting == 'full'):
            vPv = dot(vg, dot(P, vg))

        sigma0 = sqrt(vPv / (nobs - npar))

        # Standard deviations of estimated transformation parameters
        sT = sigma0 * numpy.sqrt(numpy.diag(Q))

        # Compute covariance matrix of residuals (Qv)
        Qv = None
        if (weighting == 'identity'):
            Qv = numpy.eye(nobs) - dot(dot(A, Q), A.T)
        elif (weighting == 'diagonal'):
            Qv = numpy.diag([snx.param[i].sigma**2 for i in ind1]) - dot(dot(A, Q), A.T)
        else:
            Qv = snx.Q[numpy.ix_(ind1,ind1)] - dot(dot(A, Q), A.T)

        # Compute normalized residuals
        vgn = vg / sigma0 / numpy.sqrt(numpy.diag(Qv))

        # Compute E,N,H [normalized] residuals and WRMS
        vl = numpy.zeros((nobs))
        vln = numpy.zeros((nobs))
        w = numpy.zeros(3)
        d = numpy.zeros(3)
        for i in range(nsta):
            vl[3*i+0:3*i+3] = dot(R[i], vg[3*i+0:3*i+3])
            Ql = dot(R[i], dot(Qv[3*i+0:3*i+3][:,3*i+0:3*i+3], R[i].T))
            vln[3*i+0:3*i+3] = vl[3*i+0:3*i+3] / sigma0 / numpy.sqrt(numpy.diag(Ql))
            w = w + vl[3*i+0:3*i+3]**2 / numpy.diag(Ql)
            d = d + 1 / numpy.diag(Ql)
        w = numpy.sqrt(w / d)

        # Reshape residual tables
        vg.resize(nsta, 3)
        vgn.resize(nsta, 3)
        vl.resize(nsta, 3)
        vln.resize(nsta, 3)

        # Convert raw residuals and WRMS into mm
        vg = 1000 * vg
        vl = 1000 * vl
        w  = 1000 * w

        # Write log file
        if (log != None):
            if (append):
                f = open(log, 'a')
            else:
                f = open(log, 'w')

            # Write header
            f.write('================================================================================\n')
            f.write('sinex::loadest_wrt : {0} vs {1}\n'.format(snx.file, ref.file))
            f.write('================================================================================\n')
            f.write('\n')

            # Write main options and statistics
            f.write('Number of common points : {0}\n'.format(nsta))
            f.write('Number of parameters    : {0}\n'.format(npar))
            f.write('Weighting               : {0}\n'.format(weighting))
            f.write('Variance factor         : {0:.8e}\n'.format(sigma0))
            f.write('WRMS East               : {0:8.3f} mm\n'.format(w[0]))
            f.write('WRMS North              : {0:8.3f} mm\n'.format(w[1]))
            f.write('WRMS Up                 : {0:8.3f} mm\n'.format(w[2]))
            f.write('\n')

            # Degree 1 load -> geocenter conversion factors (from P. Gegout's load Love numbers)
            kz  = (0.1285877758E+01 + 2*0.8960817937E+00)/(3*rho)*sqrt(3/(4*pi))
            kxy = kz/sqrt(2)

            # Write estimated parameters and formal errors
            f.write('Estimated parameters :\n')
            f.write('----------------------\n')
            f.write('\n')

            if (center == 'CM'):
                f.write('RX  = {0:15.8e} +/- {1:15.8e} mas\n'.format(T[0], sT[0]))
                f.write('RY  = {0:15.8e} +/- {1:15.8e} mas\n'.format(T[1], sT[1]))
                f.write('RZ  = {0:15.8e} +/- {1:15.8e} mas\n'.format(T[2], sT[2]))
                f.write('C10 = {0:15.8e} +/- {1:15.8e} kg*m^-2 <=> TZ = {2:15.8e} +/- {3:15.8e} m\n'.format(T[3], sT[3], kz *T[3], kz *sT[3]))
                f.write('C11 = {0:15.8e} +/- {1:15.8e} kg*m^-2 <=> TX = {2:15.8e} +/- {3:15.8e} m\n'.format(T[4], sT[4], kxy*T[4], kxy*sT[4]))
                f.write('S11 = {0:15.8e} +/- {1:15.8e} kg*m^-2 <=> TY = {2:15.8e} +/- {3:15.8e} m\n'.format(T[5], sT[5], kxy*T[5], kxy*sT[5]))
            elif (center in ['CF', 'CN']):
                f.write('TX  = {0:15.8e} +/- {1:15.8e} mm\n'.format(T[0], sT[0]))
                f.write('TY  = {0:15.8e} +/- {1:15.8e} mm\n'.format(T[1], sT[1]))
                f.write('TZ  = {0:15.8e} +/- {1:15.8e} mm\n'.format(T[2], sT[2]))
                f.write('RX  = {0:15.8e} +/- {1:15.8e} mas\n'.format(T[3], sT[3]))
                f.write('RY  = {0:15.8e} +/- {1:15.8e} mas\n'.format(T[4], sT[4]))
                f.write('RZ  = {0:15.8e} +/- {1:15.8e} mas\n'.format(T[5], sT[5]))
                f.write('C10 = {0:15.8e} +/- {1:15.8e} kg*m^-2 <=> TZ = {2:15.8e} +/- {3:15.8e} m\n'.format(T[6], sT[6], kz *T[6], kz *sT[6]))
                f.write('C11 = {0:15.8e} +/- {1:15.8e} kg*m^-2 <=> TX = {2:15.8e} +/- {3:15.8e} m\n'.format(T[7], sT[7], kxy*T[7], kxy*sT[7]))
                f.write('S11 = {0:15.8e} +/- {1:15.8e} kg*m^-2 <=> TY = {2:15.8e} +/- {3:15.8e} m\n'.format(T[8], sT[8], kxy*T[8], kxy*sT[8]))

            for l in range(2, lmax+1):
                if (center == 'CM'):
                    ind = (l-1)*(l+1)+3
                elif (center in ['CF', 'CN']):
                    ind = (l-1)*(l+1)+6
                f.write('C{0}0 = {1:15.8e} +/- {2:15.8e} kg*m^-2\n'.format(l, T[ind], sT[ind]))
                for m in range(1, l+1):
                    if (center == 'CM'):
                        ind = (l-1)*(l+1)+2*m+2
                    elif (center in ['CF', 'CN']):
                        ind = (l-1)*(l+1)+2*m+5
                    f.write('C{0}{1} = {2:15.8e} +/- {3:15.8e} kg*m^-2\n'.format(l, m, T[ind], sT[ind]))
                    ind = ind+1
                    f.write('S{0}{1} = {2:15.8e} +/- {3:15.8e} kg*m^-2\n'.format(l, m, T[ind], sT[ind]))
            f.write('\n')

            # Write residuals
            f.write('Residuals :\n')
            f.write('-----------\n')
            f.write('\n')
            f.write('             |  Raw residuals (mm)  | Normalized residuals |\n')
            f.write('-------------|----------------------|----------------------|\n')
            f.write('code pt soln |    E      N      H   |    E      N      H   |\n')
            f.write('-------------|----------------------|----------------------|\n')
            for i in range(nsta):
                f.write('{0} {1} {2} |'.format(code[i], pt[i], soln[i]))
                f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} |'.format(vl[i]))

                f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} |\n'.format(vln[i]))
            f.write('-------------|----------------------|----------------------|\n')
            f.write('\n')

            f.close()

        return (T, sigma0**2*Q, code, pt, soln, vl, vln, w)
