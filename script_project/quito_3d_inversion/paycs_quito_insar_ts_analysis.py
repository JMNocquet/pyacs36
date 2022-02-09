#!/usr/bin/env python

# import
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pyacs.lib.astrotime
from datetime import datetime
import numpy as np
from pyacs.gts.Gts import Gts

# path
path_InSAR='/Users/nocquet/Dropbox/raw_data/'

###############################################################################
# LOAD DATA INSAR JOHANN TIME SERIES
###############################################################################

matlab_mat=path_InSAR+'PS_StaMPS.mat'
InSAR=sio.loadmat(matlab_mat)['M']

matlab_mat2=path_InSAR+'Date_Acq.mat'
DATA_TIME=sio.loadmat(matlab_mat2)['A']

matlab_mat3=path_InSAR+'PS_Intf.mat'
DATE=sio.loadmat(matlab_mat3)['Intf_PS']

print("shape InSAR", InSAR.shape)
print("shape DATE", DATE.shape)

###############################################################################
# GET DATES IN PYACS GTS FORMAT
###############################################################################

np_ord = np.array(DATE[:,0],dtype=int)
np_dec_year = pyacs.lib.astrotime.datetime2decyear (list(map(datetime.fromordinal,np_ord))) -1.

###############################################################################
# MAIN LOOP ON PS POINTS
###############################################################################

nps = InSAR.shape[0]
#nps = 10

RESULTS = np.zeros((nps,14))

# fills longitudes and latitudes
RESULTS[:,0] = InSAR[:nps,0]
RESULTS[:,1] = InSAR[:nps,1]

for i in np.arange(nps):

    print("-- # %05d / %05d %03d %%" % (i,nps,int(i/nps*100.)) )
    
    # initialize a ts
    ts=Gts()
    ts.data =  np.zeros( (np_dec_year.shape[0],10) )
    
    ts.data[:,0] = np_dec_year
    ts.data[:,1] = InSAR[i,2:20]
    ts.data[:,2] = InSAR[i,2:20]
    ts.data[:,3] = InSAR[i,2:20]
    ts.data[:,4] = 1.
    ts.data[:,5] = 1.
    ts.data[:,6] = 1.

    # Full ts L2 norm
    dts = ts.detrend()
    vel = dts.velocity
    wrms = dts.wrms()
    RESULTS[i,2] = vel[0]
    RESULTS[i,3] = wrms[0]
    
    # L2 norm with two offsets
    dts = ts.add_offsets_dates([1994.5,1996.0]).detrend()
    vel = dts.velocity
    wrms = dts.wrms()
    RESULTS[i,4] = vel[0]
    RESULTS[i,5] = wrms[0]
    
    offset1 = dts.offsets_values[0,1]
    offset2 = dts.offsets_values[1,1]

    RESULTS[i,6] = offset1
    RESULTS[i,7] = offset2
    
    # L1 norm full ts

    dts = ts.detrend(method='L1')
    vel = dts.velocity
    wrms = dts.wrms()
    RESULTS[i,8] = vel[0]
    RESULTS[i,9] = wrms[0]
    
    # L1 norm with two offsets
    try:
        dts = ts.add_offsets_dates([1994.5,1996.0]).detrend(method='L1')
        vel = dts.velocity
        wrms = dts.wrms()
        RESULTS[i,10] = vel[0]
        RESULTS[i,11] = wrms[0]
    
        offset1 = dts.offsets_values[0,1]
        offset2 = dts.offsets_values[1,1]
    
        RESULTS[i,12] = offset1
        RESULTS[i,13] = offset2
    except:
        RESULTS[i,10] = 0.
        RESULTS[i,11] = 0.
    
        offset1 = 0.
        offset2 = 0.
    
        RESULTS[i,12] = offset1
        RESULTS[i,13] = offset2
        
    
    
print("-- Saving PYACS_LOS_VEL.dat")
np.savetxt('PYACS_LOS_VEL.dat', RESULTS, fmt="%12.7lf %12.7lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf", \
           header="long. lat. V_L2full rms_L2full V_L2 rms_L2 offset1 offset2 V_L1full rms_L1full V_L1 rms_L1 offset1 offset2")

# for i in np.arange(InSAR[:,2:].shape[1]):
#     time_in_years=7.34
#     print("-- Plotting # %d %d " % (i,DATE[i,0]))
#     fig0 = plt.figure(0,figsize=(8, 7))
#     ax0 = fig0.add_subplot(111, aspect='equal')
#     plt.ylim(-0.8,0.4)
#     plt.xlim(-79,-78.05)
#     #Background
#     #tif=TIFF.open(os.path.join(path_bkg,'geotiff/ecuador_bathy_topo_new.tif'),mode='r')
#     #topo = tif.read_image()
#     #plt.imshow(topo,cmap=plt.cm.gray,extent=[-82.1,-73.9,-5.9,3.9])
#     #InSAR
#     plt.scatter(InSAR[:,0],InSAR[:,1],c=InSAR[:,i]/time_in_years,s=0.1,cmap='jet',vmin=-3,vmax=3)
#     cbar=plt.colorbar(aspect=50)
#     cbar.set_label('Cumulated Displacement/7 years (mm/yr)', rotation=90, fontweight='bold')
#     plt.plot(-78.355,-0.221,'k*',markersize=7,label='InSAR reference')
#     plt.legend(loc='upper left',fontsize='8')
#     #Scale
#     ax0.add_patch(patches.Rectangle((-78.18, -0.77), 0.05, 0.01, facecolor='black', edgecolor='black'))
#     ax0.add_patch(patches.Rectangle((-78.13, -0.77), 0.05, 0.01, facecolor='white', edgecolor='black'))
#     ax0.text(-78.18,-0.75, '11 km', fontsize=10, weight='bold')
#     #Drawing North
#     plt.plot(range(10))
#     newax = fig0.add_axes([0.685, 0.82, 0.05, 0.05])
#     #north=plt.imread(os.path.join(path_bkg,'Narrow.png'))
#     #newax.imshow(north)
#     newax.axis('off')
#     #plt.show()
#         