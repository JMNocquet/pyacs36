#!/usr/bin/env python

# import
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
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


###############################################################################
# GET DATES IN PYACS GTS FORMAT
###############################################################################

np_ord = np.array(DATE[:,0],dtype=int)
np_dec_year = pyacs.lib.astrotime.datetime2decyear (list(map(datetime.fromordinal,np_ord))) -1.

###############################################################################
# LOAD PYACS RESULTS
###############################################################################

INSAR = np.genfromtxt(path_InSAR+'PYACS_LOS_VEL.dat')

###############################################################################
# FIGURE STAMPS RESULTS
###############################################################################

fig0 = plt.figure(0,figsize=(15, 10))
#plt.legend(loc='upper left',fontsize='8')
#Scale
#ax0.add_patch(patches.Rectangle((-78.18, -0.77), 0.05, 0.01, facecolor='black', edgecolor='black'))
#ax0.add_patch(patches.Rectangle((-78.13, -0.77), 0.05, 0.01, facecolor='white', edgecolor='black'))
#ax0.text(-78.18,-0.75, '11 km', fontsize=10, weight='bold')
#Drawing North
#plt.plot(range(10))
#newax = fig0.add_axes([0.685, 0.82, 0.05, 0.05])
#north=plt.imread(os.path.join(path_bkg,'Narrow.png'))
#newax.imshow(north)
#newax.axis('off')  


# PLOT 1 raw_ps_mean
###############################################################################

ax0 = fig0.add_subplot(2,3,1, aspect='equal',)
plt.ylim(-0.7,0.3);plt.xlim(-78.95,-78.25)
#Background
#tif=TIFF.open(os.path.join(path_bkg,'geotiff/ecuador_bathy_topo_new.tif'),mode='r')
#topo = tif.read_image()
#plt.imshow(topo,cmap=plt.cm.gray,extent=[-82.1,-73.9,-5.9,3.9])
#InSAR
plt.scatter(INSAR[:,0],INSAR[:,1],c=INSAR[:,2],s=0.1,cmap='jet',vmin=-3,vmax=3)
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar=plt.colorbar(aspect=50,cax=cax)
cbar.set_label('1993-2000 mean los rate (mm/yr)', rotation=90)
plt.plot(-78.355,-0.221,'k*',markersize=7,label='InSAR reference')

# PLOT 2 elevation
###############################################################################

ax1 = fig0.add_subplot(2,3,2, aspect='equal')
plt.ylim(-0.7,0.3);plt.xlim(-78.95,-78.25)
#Background
#tif=TIFF.open(os.path.join(path_bkg,'geotiff/ecuador_bathy_topo_new.tif'),mode='r')
#topo = tif.read_image()
#plt.imshow(topo,cmap=plt.cm.gray,extent=[-82.1,-73.9,-5.9,3.9])
#InSAR
plt.scatter(InSAR[:,0],InSAR[:,1],c=InSAR[:,-1],s=0.1,cmap='jet')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar=plt.colorbar(aspect=50,cax=cax)
cbar.set_label('height (m)', rotation=90)

# PLOT 3 uncertainties
###############################################################################

ax2 = fig0.add_subplot(2,3,3, aspect='equal')
plt.ylim(-0.7,0.3);plt.xlim(-78.95,-78.25)
#Background
#tif=TIFF.open(os.path.join(path_bkg,'geotiff/ecuador_bathy_topo_new.tif'),mode='r')
#topo = tif.read_image()
#plt.imshow(topo,cmap=plt.cm.gray,extent=[-82.1,-73.9,-5.9,3.9])
#InSAR
plt.scatter(InSAR[:,0],InSAR[:,1],c=InSAR[:,-2],s=0.1,cmap='jet')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar=plt.colorbar(aspect=50,cax=cax)
cbar.set_label('LOS rate uncertainty (mm/yr)', rotation=90)

# PLOT 4 corrected LOS
###############################################################################

ax3 = fig0.add_subplot(2,3,5, aspect='equal')
plt.ylim(-0.7,0.3);plt.xlim(-78.95,-78.25)
#Background
#tif=TIFF.open(os.path.join(path_bkg,'geotiff/ecuador_bathy_topo_new.tif'),mode='r')
#topo = tif.read_image()
#plt.imshow(topo,cmap=plt.cm.gray,extent=[-82.1,-73.9,-5.9,3.9])
#InSAR

lindex = np.where(InSAR[:,-2]<=1.5)
plt.scatter(INSAR[lindex,0],INSAR[lindex,1],c=INSAR[lindex,4],s=0.2,alpha=0.1,cmap='jet',  vmin=-3.,vmax=3.)
lindex = np.where(InSAR[:,-2]<=.9)
plt.scatter(INSAR[lindex,0],INSAR[lindex,1],c=INSAR[lindex,4],s=0.2,cmap='jet',  vmin=-3.,vmax=3.)
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar=plt.colorbar(aspect=50,cax=cax)
cbar.set_label('1996-2000 LOS rate (mm/yr)', rotation=90)

# PLOT 5 transient
###############################################################################

ax4 = fig0.add_subplot(2,3,4, aspect='equal')
plt.ylim(-0.7,0.3);plt.xlim(-78.95,-78.25)
#Background
#tif=TIFF.open(os.path.join(path_bkg,'geotiff/ecuador_bathy_topo_new.tif'),mode='r')
#topo = tif.read_image()
#plt.imshow(topo,cmap=plt.cm.gray,extent=[-82.1,-73.9,-5.9,3.9])
#InSAR
lindex = np.where(InSAR[:,-2]<=.9)
plt.scatter(INSAR[:,0],INSAR[:,1],c=INSAR[:,6]+INSAR[:,7],s=0.1,cmap='jet', vmin=-25.,vmax=25.)
#plt.scatter(INSAR[lindex,0],INSAR[lindex,1],c=INSAR[lindex,6]+INSAR[lindex,7],s=0.2,cmap='jet',vmin=-15.,vmax=15.)
divider = make_axes_locatable(ax4)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar=plt.colorbar(aspect=50,cax=cax)
cbar.set_label('LOS 1995/4/11-1996/10/24 (mm)', rotation=90)

INSAR_TRANSIENT = np.zeros( (INSAR.shape[0] , 3 ) )
INSAR_TRANSIENT[:,:2] = INSAR[:,:2]
INSAR_TRANSIENT[:,2] = INSAR[:,6]+INSAR[:,7]

lindex = np.where(InSAR[:,-2]<=.9)

np.savetxt('INSAR_TRANSIENT_1995_1996.dat', INSAR_TRANSIENT[lindex], fmt = "%12.7lf %12.7lf %5.4lf  ")

 

###############################################################################
# PLOT 5 extract transient

lon_ref1 = -78.39
lat_ref1 = -0.37
dc1 = 0.15
lindex_silalo1 = np.where( np.sqrt((INSAR[:,0]-lon_ref1)**2 + (INSAR[:,1]-lat_ref1)**2) < dc1 )

lon_ref2 = -78.4
lat_ref2 = -0.16
dc2 = 0.08
lindex_silalo2 = np.where( np.sqrt((INSAR[:,0]-lon_ref2)**2 + (INSAR[:,1]-lat_ref2)**2) < dc2 )

lindex_silalo = ( np.array(list(set( lindex_silalo1[0].tolist() + lindex_silalo2[0].tolist()))) ,)

ax5 = fig0.add_subplot(2,3,6, aspect='equal')
plt.ylim(-0.7,0.3);plt.xlim(-78.95,-78.25)

#plt.scatter(INSAR[lindex_silalo,0],INSAR[lindex_silalo,1],c=INSAR[lindex_silalo,6]+INSAR[lindex_silalo,7],s=0.2,cmap='jet',vmin=-15.,vmax=15.)

INSAR_FINAL = np.copy(INSAR[:,:4])
INSAR_FINAL[:,3] = InSAR[:,-2]
INSAR_FINAL[:,2] = (INSAR_FINAL[:,2] + INSAR[:,4] )/2.
 
INSAR_FINAL[lindex_silalo,2] = INSAR[lindex_silalo,4]


lindex = np.where(InSAR[:,-2]<=1.5)
plt.scatter(INSAR_FINAL[lindex,0],INSAR_FINAL[lindex,1],c=INSAR_FINAL[lindex,2],alpha=0.1,s=0.1,cmap='jet', vmin=-3.,vmax=3.)

lindex = np.where(InSAR[:,-2]<=.9)
plt.scatter(INSAR_FINAL[lindex,0],INSAR_FINAL[lindex,1],c=INSAR_FINAL[lindex,2],s=0.2,cmap='jet',  vmin=-3.,vmax=3.)

divider = make_axes_locatable(ax5)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar=plt.colorbar(aspect=50,cax=cax)
cbar.set_label('1993-2000 corrected LOS rate (mm/yr)', rotation=90)


plt.subplots_adjust(wspace = 0.5)
#plt.tight_layout()
plt.show()

np.savetxt('INSAR_FINAL.dat', INSAR_FINAL[lindex], fmt = "%12.7lf %12.7lf %5.4lf %5.4lf  ")
