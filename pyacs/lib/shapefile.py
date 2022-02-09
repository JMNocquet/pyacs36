"""
Various subroutine to convert pyacs results to shapefiles
"""

###############################################################################
def static_slip_to_shapefile(static_slip_file,shp_name,dis_type='rec',verbose=False):
###############################################################################
    """
    Converts pyacs/pyeq slip solution into a shapefile
    
    :param static_slip_file: output from pyeq_static_inversion.py (*sol_slip.dat)
    :param shp_name: output shapefile name
    :param dis_type: either 'rec' (rectangular dislocation) or 'tde' (triangular dislocation element)
    :param verbose: verbose mode
    
    """
    
    import shapefile
    from pyeq.lib import eq_disloc_3d as DL
    import numpy as np


    lfaults=[]
    lrecord=[]

    # rectangular or triangular dislocation
    if dis_type == 'tde':
        if verbose:
            print("-- Triangular dislocations output")
        TRIANGLE=True
    else:
        if verbose:
            print("-- Rectangular dislocations output")
        TRIANGLE=False

    # reads file
    try:
        
        SLIP = np.array( np.asmatrix( np.genfromtxt(static_slip_file,comments='#') ) )
        
    except:
        raise IOError("!!! Could not read ",static_slip_file)

    print('SLIP.shape' , SLIP.shape )
    
    # loop on fault elements
    
    for i in np.arange(SLIP.shape[0]):

        (rdis_long, rdis_lat, rdis_depth, rdis_length, rdis_width,  rdis_area, _ratio_rdis_tdis,strike,dip,\
         _centroid_long, _centroid_lat, centroid_depth,\
         tdis_long1,tdis_lat1,_tdis_depth1,tdis_long2,tdis_lat2,_tdis_depth2,tdis_long3,tdis_lat3,_tdis_depth3,_tdis_area,\
         _rake_1, _slip_1 ,_rake_2 ,  _slip2  ,  slip) \
         =SLIP[i,:]

    
        if TRIANGLE:
            lfaults.append([ [tdis_long1,tdis_lat1], [tdis_long2,tdis_lat2], [tdis_long3,tdis_lat3] ])
            lrecord.append([i,centroid_depth,slip,0])
    
        else:
        
            # creates a dislocation object
            disloc=DL.Dislocation(i,rdis_long,rdis_lat,rdis_depth, strike, dip, rdis_length, rdis_width,rdis_area, 0, 0)
            # get the corners
            (X1,X2,X3,X4)=disloc.corners(coor_type='geo')
            lfaults.append([ [X1[0],X1[1]], [X2[0],X2[1]], [X3[0],X3[1]], [X4[0],X4[1]], [X1[0],X1[1]] ])
            lrecord.append([i,rdis_depth,slip,0])

    print("-- ",len(lfaults)," polygons read")

    ###################################################################
    # WRITES SHAPEFILES
    ###################################################################
    
    ###################################################################
    # INITIALIZE SHAPEFILE
    ###################################################################
    
    # Make a polygon shapefile

    w = shapefile.Writer(shp_name,shapeType=shapefile.POLYGON)
    w.field('ID','I','40')
    w.field('i_subfault','F','40')
    w.field('depth_top_disloc','F','40')
    w.field('slip','F','40')
    w.field('rake','F','40')
    
    
    ###################################################################
    # LOOP ON FAULTS
    ###################################################################
    
    for i in np.arange(len(lfaults)):
        fault=lfaults[i]
        w.poly([fault])
        [index,depth,slip,rake]=lrecord[i]
        w.record(str(i),index,depth,slip,rake)
        i=i+1
    
    ###################################################################
    # SAVE SHAPEFILE
    ###################################################################
    
    shapefile=shp_name+'.shp'
    print(("-- saving shapefile %s " % (shapefile)))
    w.close()
    
    ###################################################################
    # SAVE .PRJ FILE
    ###################################################################
    
    prj = open("%s.prj" % shp_name, "w") 
    epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]' 
    prj.write(epsg) 
    prj.close()
    
    ###################################################################
    # WRITES GMT PSXY FILES
    # This file can be then plotted with
    #  psxy lima_simple_sigma_010_dc_050_m0_000_sol_coupling.gmt -V -R-81/-72/-16/-7 -JM14 -L -m -W0.2/0  -Clut_coupling.cpt > test.ps
    # triangles have not been tested yet
    ###################################################################

    gmtfile=shp_name+'.gmt'

    print(("-- saving gmt file %s " % gmtfile.split('/')[-1]))
    f=open(gmtfile,'w')
    
    for i in np.arange(len(lfaults)):
        [index,depth,slip,rake]=lrecord[i]
        f.write('> -Z%.3lf\n'%slip)
        fault=lfaults[i]
        for xy in fault:
            f.write("%10.3lf %10.3lf\n" %(xy[0],xy[1]))
    
    f.write('>\n')
    f.close()
    

###############################################################################
def psvelo_to_shapefile(psvelo_file,shp_name,verbose=False):
###############################################################################
    """
    Converts a psvelo GMT file into shapefile
    
    :param psvelo: GMT psvelofile
    :param shp_name: output shapefile name
    :param verbose: verbose mode
    
    """

    import shapefile
    from pyacs.lib.vel_field import Velocity_Field as VF
    import numpy as np

    ###################################################################
    # READS GMT FILE
    ###################################################################
    vel=VF.read(file_name=psvelo_file,verbose=verbose)
    
    # case empty file
    if vel.nsites() == 0:
        from colors import red
        print( red("[PYACS WARNING]: empty psvelo file: %s " % psvelo_file ))
        return
    
    ###################################################################
    # INITIALIZE SHAPEFILE
    ###################################################################
    
    # Make a point shapefile
    w = shapefile.Writer(shp_name,shapeType=shapefile.POINT)
    w.field('longitude','C','40')
    w.field('latitude','C','40')
    w.field('Ve','C','40')
    w.field('Vn','C','40')
    w.field('S_Ve','C','40')
    w.field('S_Vn','C','40')
    w.field('S_Ven','C','40')
    w.field('name','C','40')
    w.field('name_velocity','C','40')
    
    
    ###################################################################
    # LOOP ON SITES
    ###################################################################
    
    lsite=vel.lcode()
    
    for site in lsite:
        
        if verbose:
            print('-- adding ',site,' to shapefile')
        
        M=vel.site(site)
        w.point(M.lon,M.lat)
        name_velocity=("%s %5.2lf" % (M.code,np.sqrt(M.Ve**2+M.Vn**2)))
        w.record(M.lon,M.lat,M.Ve,M.Vn,M.SVe,M.SVn,'0',M.code,name_velocity)
    
    
    ###################################################################
    # SAVE SHAPEFILE
    ###################################################################
    
    shapefile=shp_name+'.shp'
    w.close()
    if verbose:
        print(("-- saving shapefile %s " % (shapefile)))

    ###################################################################
    # SAVE .PRJ FILE
    ###################################################################
    
    prj = open("%s.prj" % shp_name, "w") 
    epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]' 
    prj.write(epsg) 
    prj.close()


###############################################################################
def pyeblock_fault(pyeblock_file,shp_name, verbose=False):
###############################################################################
    """
    Converts a pyeblock fault file into a shapefile
    
    :param pyeblock_fault: pyeblock fault file
    :param shp_name: output shapefile name
    :param verbose: verbose mode
    
    """
    
    import shapefile
    import numpy as np


    lfaults=[]
    lrecord=[]

    # reads a pyeblock format
    try:
        FAULTS      = np.array( np.asmatrix( np.genfromtxt(pyeblock_file, usecols=range(8) , comments='#') ) )
        FAULTS_POLY = np.array(np.asmatrix( np.genfromtxt(pyeblock_file,usecols=(8,9), dtype=str ) ))
        
    except:
        raise IOError("!!! Could not read ",pyeblock_file)

    
    # loop on fault elements
    
    ###################################################################
    # WRITES SHAPEFILES
    ###################################################################
    
    ###################################################################
    # INITIALIZE SHAPEFILE
    ###################################################################
    
    # Make a polyline shapefile

    w = shapefile.Writer(shp_name,shapeType=shapefile.POLYLINE)
    w.field('ID','I','40')
    w.field('top','F','40')
    w.field('bottom','F','40')
    w.field('dip','F','40')
    w.field('left_block','C','40')
    w.field('right_block','C','40')
    
    
    ###################################################################
    # LOOP ON FAULTS
    ###################################################################
    
    for i in np.arange(FAULTS.shape[0]):
        print( FAULTS[i,1:3] )
        w.line( [ [FAULTS[i,1:3] , FAULTS[i,3:5]] ] )
        w.record(str(i), FAULTS[i,5],FAULTS[i,6],FAULTS[i,7],FAULTS_POLY[i,0],FAULTS_POLY[i,1] )
    
    ###################################################################
    # SAVE SHAPEFILE
    ###################################################################
    
    shapefile=shp_name+'.shp'
    print(("-- saving shapefile %s " % (shapefile)))
    w.close()
    
    ###################################################################
    # SAVE .PRJ FILE
    ###################################################################
    
    prj = open("%s.prj" % shp_name, "w") 
    epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]' 
    prj.write(epsg) 
    prj.close()
    
    ###################################################################
    # WRITES GMT PSXY FILES
    # This file can be then plotted with
    #  psxy lima_simple_sigma_010_dc_050_m0_000_sol_coupling.gmt -V -R-81/-72/-16/-7 -JM14 -L -m -W0.2/0  -Clut_coupling.cpt > test.ps
    # triangles have not been tested yet
    ###################################################################

    gmtfile=shp_name+'.gmt'

    print(("-- saving gmt file %s " % gmtfile.split('/')[-1]))
    f=open(gmtfile,'w')
    
    for i in np.arange(FAULTS.shape[0]):
        f.write('> -Z%d\n'% i)
        f.write("%10.3lf %10.3lf\n" %( FAULTS[i,1] , FAULTS[i,2] ))
        f.write("%10.3lf %10.3lf\n" %( FAULTS[i,3] , FAULTS[i,4] ))
    
    f.write('>\n')
    f.close()



