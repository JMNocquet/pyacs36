#!/usr/bin/env python

import os
from pyacs.lib import syscmd
import numpy as np
from pyeq.lib.eq_disloc_3d import Dislocation

# files

hgps =  '/Users/nocquet/projets/2018/inversion_quito_3d/pyeq/data_corrected_subduction/hgps_corr_subd_quito.dat'
insar = '/Users/nocquet/projets/2018/inversion_quito_3d/pyeq/data_corrected_subduction/insar_corr_subd_quito_new_500m.dat'

# name experiment

expt = 'quito_01'

GEOMETRY = True
GREEN    = True
MODEL    = True
PLOT     = True

# Make a geometry for Quito
# The geometry has two parts: 
#    - a flat decollement to explain the horizontal
#    - a ramp corresponding to the Quito inverse fault system

# First create the working directories

try:
    os.mkdir('geometry')
except:
    pass

try:
    os.mkdir('green')
except:
    pass

try:
    os.mkdir('model')
except:
    pass

###############################################################################
if GEOMETRY:
###############################################################################

    # general parameters
    
    for dip in [50]:
            
#            dip   = 30.
 
        for depth_decollement in [7,10,13]: 
        #depth_decollement = 13.
        
            width = depth_decollement / np.sin( np.radians(dip) )
            
            print("-- Geometry dip=%.1lf depth=%.1lf width=%.1lf " % (dip,depth_decollement,width))
            
            # ramp
            
            lon0 = -78.3
            lat0 =   0.20
            depth = 0.
            strike = 200.
            length = 110.
            size_subfault = 2 
        
            nstrike = int( length / size_subfault )
            ndip = int( width / size_subfault )
        
            cmd = 'cd geometry ; pyeq_make_rectangular_fault.py '\
                  + ' -corner /'+str(lon0)+'/'+str(lat0)+'/'+str(depth)\
                  + ' -length '+str(length)\
                  + ' -width '+ str(width)\
                  + ' -strike '+ str(strike)\
                  + ' -dip '+str(dip)\
                  + ' -e ' + expt\
                  + ' --nstrike '+ str(nstrike)\
                  + ' --ndip '+str(ndip)\
                  + ' --verbose'
            
            print('-- running: ',cmd)
            syscmd.getstatusoutput( cmd , verbose=False )
        
            np_ramp = np.load('geometry/'+expt+'_geometry.npy')
        
            # get info for decollement
            
            from pyeq.lib.eq_disloc_3d import Dislocation
            
            my_ramp = Dislocation( 0 , lon0, lat0, depth, strike, dip, length, width,0, 0, 0)
        
            [lon_decollement,lat_decollement,depth_decollement] = my_ramp.corners(coor_type='geo')[3]
            print([lon_decollement,lat_decollement,depth_decollement])
            # decollement
            
            width = length
            dip_decollement = 0.1
        
            size_subfault = 1000
            nstrike = int( length / size_subfault )
            ndip = int( width / size_subfault )
            
            cmd = 'cd geometry ; pyeq_make_rectangular_fault.py '\
                  + ' -corner /'+str(lon_decollement)+'/'+str(lat_decollement)+'/'+str(depth_decollement)\
                  + ' -length '+str(length)\
                  + ' -width '+ str(width)\
                  + ' -strike '+ str(strike)\
                  + ' -dip '+str(dip_decollement)\
                  + ' -e ' + expt\
                  + ' --nstrike '+ str(nstrike)\
                  + ' --ndip '+str(ndip)\
                  + ' --verbose'
            
            print('-- running: ',cmd)
            syscmd.getstatusoutput( cmd , verbose=False )
        
            np_decollement = np.load('geometry/'+expt+'_geometry.npy')
        
            # merge geometry
            
            print('-- merging geometry')
            
            np_geometry = np.insert(np_ramp,0,np_decollement,axis=0)
            
            print('-- geometry ', np_geometry.shape)
            # write geometry
            
            names=['rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
                   'rdis_area','ratio_rdis_tdis','strike','dip',\
                   'centroid_long','centroid_lat','centroid_depth',\
                   'tdis_long1','tdis_lat1','tdis_depth1',\
                   'tdis_long2','tdis_lat2','tdis_depth2',\
                   'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']
             
            np.save('geometry/'+expt+'_geometry.npy', np_geometry)
             
            GEOMETRY_TXT=np.zeros((np_geometry.shape[0],np_geometry.shape[1]+1))
            GEOMETRY_TXT[:,0]=np.arange(np_geometry.shape[0])
            GEOMETRY_TXT[:,1:] = np_geometry
             
            Hdtype={'names':['    rdis_long','   rdis_lat','rdis_depth',\
                            'rdis_length','rdis_width',' rdis_area','ratio_rdis_tdis',\
                            '   strike','       dip',\
                            'centroid_long','centroid_lat','centroid_depth',\
                            'tdis_long1', ' tdis_lat1','tdis_depth1', \
                            'tdis_long2',' tdis_lat2','tdis_depth2', \
                            'tdis_long3',' tdis_lat3','tdis_depth3', \
                            ' tdis_area'],\
                     'formats':['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']}
         
            header_cols=' '.join(Hdtype['names'])
            format_header='%04d %10.5lf %10.5lf      %6.2lf \
                 %6.2lf     %6.2lf    %6.2lf          %6.2lf\
                %6.2lf %10.2lf \
               %10.5lf   %10.5lf         %6.2lf \
            %10.5lf %10.5lf      %6.2lf \
            %10.5lf %10.5lf      %6.2lf \
            %10.5lf %10.5lf      %6.2lf \
               %6.2lf '
             
            np.savetxt('geometry/'+expt+'_geometry.dat', GEOMETRY_TXT, fmt=format_header, header=header_cols)
        
            
        
            ###############################################################################
            if GREEN:
            ###############################################################################
                # make the Green's function
                
                
                cmd = 'cd green ; pyeq_make_green.py '\
                      + ' -g ../geometry/'+expt+'_geometry.npy'\
                      + ' -gps_h '+hgps\
                      + ' -insar '+insar\
                      + ' -e '+expt\
                      + ' --verbose '
                
                print('-- running: ',cmd)
                from time import time
                t0 = time()
            
                syscmd.getstatusoutput( cmd , verbose=False )
            
                print('time = ', time() - t0 )
            
            ###############################################################################
            if MODEL:   
            ###############################################################################
                
                # run the model
                
#                for vh in [3,4.5,6,8]:
                for vh in [ 8. ]:
                    for mm0 in [0.,1.]:
                        for sigma in [1.,10.]:
                            for dc in [2.,4.]:
            #                    print('os.path ' , os.getcwd())
            #            vh = 4.5
            #            dc = 4.
            #            sigma = 1.
                                rake_type = 'constant'
                                rake_constraint = 1.
                                rake_value = 106.
                                # s_ramp = vh / cos(dip)
                                max_slip = ( vh + 0.1 ) / np.cos( np.radians( dip) )
                                # s_ramp = Vh
                                max_slip = vh
                                m0 = mm0 * max_slip
                                wh = 1.
                                wi = 1.
                                insar_opt = ' --insar '
                                str_vh = ("0/%.1lf/%.1lf" % (vh-0.1,vh+0.1) )
                            #    idx_decollement = np.genfromtxt('geometry/'+expt+'_geometry.dat').shape[0]-1
                                idx_decollement = 0
                                print('-- idx_decollement: ',idx_decollement)
                            
                                # outdir
                                
                                outdir = ("dip_%2d_depth_%02d_vh_%02d_m0_%1d_sigma_%02d_dc_%1d" % (dip,int(depth_decollement),vh,int(m0/max_slip),int(sigma),int(dc)) )
                            
                                cmd = 'cd model ; pyeq_static_inversion.py '\
                                      + ' -input_npz ../green/'+expt+'_input.npz '\
                                      + ' -dc '+str(dc)\
                                      + ' -sigma '+str(sigma)\
                                      + ' -m0 '+str(m0)\
                                      + ' -rake_type '+ rake_type\
                                      + ' -rake_value ' + str(rake_value)\
                                      + ' -rake_constraint '+ str(rake_constraint)\
                                      + ' -max_slip '+ str(max_slip)\
                                      + ' -wh '+str(wh) \
                                      + ' -wi '+str(wi) \
                                      + ' '+insar_opt \
                                      + ' -e '+outdir\
                                      + ' --verbose'\
                                      + ' --c '+str_vh\
                                      + ' ; cd ..'
                                     
                                print('-- running: ',cmd)
                                print( syscmd.getstatusoutput( cmd , verbose=True ) )
                                
                                ###############################################################################
                                if PLOT:
                                ###############################################################################
                                
                                    from glob import glob
                                    verbose=False
                                    
                                    # find the latest directory
                                    wdir = max(glob(os.path.join('model', '*/')), key=os.path.getmtime)
                                
                                    print('-- making plot for dir: ',wdir)
                                    
                                    cmd = 'cd '+wdir
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    os.chdir(wdir)
                                    
                                    # read results
                                
                                    sum_dat = glob( '*sum.dat' )[0]
                                
                                    import re
                                    
                                    with open(sum_dat) as f:
                                        content = f.readlines()
                                        # you may also want to remove whitespace characters like `\n` at the end of each line
                                        content = [x.strip() for x in content] 
                                
                                    var_red_gpsh  = list(filter(lambda x: re.search(r'variance reduction hor. GPS', x), content))[0].split(':')[-1].lstrip()
                                    var_red_insar = list(filter(lambda x: re.search(r'variance reduction InSAR', x), content))[0].split(':')[-1].lstrip()
                                
                                    wrms_gpsh   = list(filter(lambda x: re.search(r'wrms h_gps', x), content))[0].split(':')[-1].lstrip()
                                    wrms_insar  = list(filter(lambda x: re.search(r'wrms InSAR', x), content))[0].split(':')[-1].lstrip()
                                    print( list(filter(lambda x: re.search(r'slip max', x), content))[0].split(':')[-1].split()[0] )
                                    slip_max  = str( list(filter(lambda x: re.search(r'slip max', x), content))[0].split(':')[-1].split()[0] )
                                    
                                    
                                    # make  plot of results using GMT
                                    
                                    
                                    # intialization
                                    try:
                                        os.remove('plot_1.ps')
                                    except:
                                        pass
                                    
                                    ps = ' >> plot_1.ps'
                                    topo = '/usr/local/geodesy/maps/geotiff/ecuador/GMRTv3_2_20161209topo.grd'
                                    bounds=' -R-78.9/-78.1/-0.7/0.3 '
                                    proj = ' -JM6 '
                                    pl = ''
                                    bm = ' -Ba.2f.1 '
                                    basemap = bounds+proj+pl
                                    lut_topo = '/usr/local/geodesy/maps/templates/ecuador_sea.cpt'
                                    cmd = 'gmt set FONT_ANNOT_PRIMARY 4p,Helvetica MAP_GRID_CROSS_SIZE_PRIMARY 0.1i MAP_ANNOT_OFFSET_PRIMARY 1p'+ps
                                    syscmd.getstatusoutput( cmd , verbose=True )
                                    
                                    cmd = ' gmt set MAP_FRAME_TYPE plain'+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    cmd = ' gmt set FORMAT_GEO_MAP D'+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = ' gmt set FONT_TITLE  8p,Helvetica,black'+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = ' gmt set MAP_TITLE_OFFSET  5p'+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt grdgradient '+topo+'  -Ggradient_topo.grd -A10/60 -Nt'
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt makecpt -Chot -N -I -T0/7/1 > slip.cpt'
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt makecpt -Cjet  -N -T-3./3./.5 > insar.cpt'
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    # plot 1: slip model
                                
                                    model_dat = glob( '*sol_slip.dat' )[0]
                                    
                                    cmd = 'pyeq_static_model_to_qgis.py -dat '+model_dat
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    model_gmt = model_dat.split('.')[0]+'.gmt'
                                    
                                    np_sol_slip = np.array( np.asmatrix( np.genfromtxt(model_dat) ))
                                    
                                    str_vh = ("%.1lf " % np_sol_slip[0,-1] )
                                    str_depth = ("%.1lf " % depth )
                                    
                                    
                                    ok=' -K '
                                    title = ' -B+t"fault slip vh='+str_vh+'mm/yr depth='+str_depth+'km slip max '+slip_max+'"'
                                    cmd = 'gmt pscoast -X1. '+title+basemap+bm+' -Dh -N1 -W1 '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    ok=' -O -K '
                                    cmd = 'gmt grdimage '+topo+' -Igradient_topo.grd -M  -C'+lut_topo+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    
                                    cmd = 'gmt psxy '+model_gmt+basemap+' -L  -W0.2/0  -Cslip.cpt'+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    cmd = 'gmt psscale -Dg-78.9/-0.8+w6/.2h -Cslip.cpt '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                 
                                    # plot 2: horizontal GPS obs & model
                                
                                    title = ' -B+t"model&obs GPS var_red='+var_red_gpsh+' wrms='+wrms_gpsh+'"'
                                    
                                    gps_obs = glob( '*_obs.dat' )[0]
                                    gps_model = glob( '*_model.dat' )[0]
                                    gps_res = glob( '*_residuals.dat' )[0]
                                
                                    cmd = 'gmt pscoast -X7.1 '+title+basemap+bm+' -Dh -N1 -W1 '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt grdimage '+topo+' -Igradient_topo.grd -M  -C'+lut_topo+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt psxy '+model_gmt+basemap+' -W0.2/0  '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                
                                    cmd = 'gmt psvelo '+gps_obs+' -G0/0/255 -A0.02/0.15/0.05 -W0,.5 -Se0.15/0.95/0 '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt psvelo '+gps_model+' -G255/0/0 -A0.02/0.15/0.05 -W.5/0/0/0 -Se0.15/0.95/0 '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt psvelo '+gps_res+' -G255/255/0 -A0.02/0.15/0.05 -W.5/0/0/0 -Se0.15/0.95/0 '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    # plot 3: InSAR model
                                
                                #    title = ' -B+t"InSAR LOS predicted var_red='+var_red_insar+'"'
                                    title = ' -B+t"InSAR LOS predicted var_red='+var_red_insar+' wrms='+wrms_insar+'"'
                                    
                                    insar_obs = glob( '*_model_insar.dat' )[0]
                                    
                                    cmd = 'gmt pscoast -X7.1 '+title+basemap+bm+' -Dh -N1 -W1 '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt grdimage '+topo+' -Igradient_topo.grd -M  -C'+lut_topo+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    cmd = 'gmt psxy '+insar_obs+basemap+' -Sc.03  -Cinsar.cpt'+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt psscale -Dg-78.9/-0.8+w6/.2h -Cinsar.cpt '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    # plot 4: InSAR obs
                                
                                    title = ' -B+t"InSAR LOS observed" '
                                    
                                    insar_obs = insar
                                    
                                    np_insar_obs = np.genfromtxt(insar_obs,usecols=(3,4,5))
                                    np_insar_obs[:,2] = np_insar_obs[:,2] * 1000. 
                                    
                                    np.savetxt('insar_obs.dat',  np_insar_obs)
                                    
                                    cmd = 'gmt pscoast -X7.1 '+title+basemap+bm+' -Dh -N1 -W1 '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt grdimage '+topo+' -Igradient_topo.grd -M  -C'+lut_topo+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    cmd = 'gmt psxy insar_obs.dat '+basemap+' -Sc.03  -Cinsar.cpt'+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt psscale -Dg-78.9/-0.8+w6/.2h -Cinsar.cpt '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    # plot 5: InSAR res
                                
                                    title = ' -B+t"InSAR LOS residuals" '
                                    
                                    insar_res = glob( '*_residuals_insar.dat' )[0]
                                    
                                    np_insar_res = np.genfromtxt(insar_res,usecols=(0,1,2))
                                    np_insar_res[:,2] = np_insar_res[:,2] 
                                    
                                    np.savetxt('insar_res.dat',  np_insar_res)
                                    
                                    cmd = 'gmt pscoast -X-21.3 -Y9.5 '+title+basemap+bm+' -Dh -N1 -W1 '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt grdimage '+topo+' -Igradient_topo.grd -M  -C'+lut_topo+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    cmd = 'gmt psxy insar_res.dat '+basemap+' -Sc.03  -Cinsar.cpt'+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt psscale -Dg-78.9/-0.8+w6/.2h -Cinsar.cpt '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    # subplot 5: model prediction over a grid
                                
                                
                                    # first generate a grid of point where model prediction will be made
                                    nx = 5
                                    ny = 2 * nx
                                    
                                    x = np.linspace(-78.9,-78.1, nx)
                                    y = np.linspace(-0.7,0.3, ny)
                                    
                                    mesh_grid = np.meshgrid(x,y)
                                    grid = np.vstack((mesh_grid[0].T.flatten(),mesh_grid[1].T.flatten())).T
                                
                                    np.savetxt('void.gmt',grid,fmt="%.8lf %.8lf")
                                    
                                #    cmd = 'tail -1 '+model_dat+' > slip_decollement_only.dat'
                                #    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    cmd = 'pyeq_model_to_disp.py -i void.gmt -model '+model_dat+' -o model_prediction_grid.dat'
                                    print("-- running: %s" % cmd)
                                    
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                
                                    m_grid = np.genfromtxt('model_prediction_grid.dat')
                                
                                
                                    # subplot horizontal displacement
                                    
                                    vh_grid = np.zeros((m_grid.shape[0],8))
                                    vh_grid[:,:4] = m_grid [:,:4]
                                    np.savetxt('vh_grid.dat',vh_grid,fmt='%.4lf %.4lf %.4lf %.4lf %.1lf %.1lf %.1lf %.1lf')
                                    
                                    title = ' -B+t"model Vh" '
                                
                                    cmd = 'gmt pscoast -X7.1 '+title+basemap+bm+' -Dh -N1 -W1 '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt psvelo '+'vh_grid.dat'+' -G0/0/255 -A0.02/0.15/0.05 -W0,.5 -Se0.1/0.95/0 '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    # subplot vertical displacement
                                
                                    nx = 5
                                    ny = 2 * nx
                                    
                                    x = np.linspace(-78.9,-78.1, nx)
                                    y = np.linspace(-0.7,0.3, ny)
                                    
                                    mesh_grid = np.meshgrid(x,y)
                                    grid = np.vstack((mesh_grid[0].T.flatten(),mesh_grid[1].T.flatten())).T
                                
                                    np.savetxt('void.gmt',grid,fmt="%.8lf %.8lf")
                                    
                                    cmd = 'pyeq_model_to_disp.py -i void.gmt -model '+model_dat+' -o model_prediction_grid.dat'
                                    print("-- running: %s" % cmd)
                                    
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                
                                    m_grid = np.genfromtxt('model_prediction_grid.dat')
                                    
                                    vu_grid = np.zeros((m_grid.shape[0],3))
                                    vu_grid[:,:2] = m_grid [:,:2]
                                    vu_grid[:,2] = m_grid [:,4]
                                    np.savetxt('vu_grid.dat',vu_grid,fmt='%.4lf %.4lf %.4lf ')
                                    
                                    title = ' -B+t"model Vu" '
                                
                                    cmd = 'gmt pscoast -X7.1 '+title+basemap+bm+' -Dh -N1 -W1 '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt grdimage '+topo+' -Igradient_topo.grd -M  -C'+lut_topo+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    cmd = 'gmt psxy '+'vu_grid.dat'+basemap+' -Sc.08  -Cinsar.cpt'+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt psscale -Dg-78.9/-0.8+w6/.2h -Cinsar.cpt '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    # locking
                                
                                    cmd = 'gmt makecpt -Chot -I -N -T0/100/10 > coupling.cpt'
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    title = ' -B+t"coupling" '
                                
                                    np_coupling = np.genfromtxt(model_dat)
                                    np_coupling[:,-1] = 100. - np_coupling[:,-1] / ( (np_sol_slip[0,-1]+0.2) / np.cos( np.radians( np_coupling[:,8] ) )) *100.  
                                    np.savetxt('coupling.dat', np_coupling)
                                
                                    cmd = 'pyeq_static_model_to_qgis.py -dat coupling.dat'
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    cmd = 'gmt pscoast -X7.1 '+title+basemap+bm+' -Dh -N1 -W1 '+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                    cmd = 'gmt grdimage '+topo+' -Igradient_topo.grd -M  -C'+lut_topo+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    cmd = 'gmt psxy coupling.gmt '+basemap+' -L  -W0.2/0  -Ccoupling.cpt'+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                
                                
                                    ok = ' -O '
                                    cmd = 'gmt psscale -Dg-78.9/-0.8+w6/.2h -Ccoupling.cpt '+basemap+ok+ps
                                    syscmd.getstatusoutput( cmd , verbose=verbose )
                                    
                                    
                                
                                    print("-- created %s/plot.ps" % (wdir))
                    
                                    os.chdir('../..')
                    
                                    
                                    #cmd = 'open -a Preview plot.ps'
                                    #syscmd.getstatusoutput( cmd , verbose=verbose )
                            