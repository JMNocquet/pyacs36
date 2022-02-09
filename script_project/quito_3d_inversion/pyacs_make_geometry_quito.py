
###############################################################################
# DECOLLEMENT - SINGLE RAMP GEOMETRY
###############################################################################

    # First make the ramp geometry
     
#     lon0 = -78.31
#     lat0 = 0.194
#     depth = 15.
#     strike = 200.
#     dip = .5
#     length = 110.
#     width  =  100.
#     size_subfault = 1000 
 
 
#     lon0 = -78.4
#     lat0 = -0.1
#     depth = 0.
#     strike = 200.
#     dip = 20.
#     length = 110.
#     width  =  25.
#     size_subfault = 2 
 
 
    nstrike = int( length / size_subfault )
    ndip = int( width / size_subfault )
 
    # test a single decollement only
 
#     lon0 = -78.46295
#     lat0 = 0.194
#     depth = 1.14
#     strike = 195.
#     dip = .1
#     length = 90.
#     width  =  100.
#     size_subfault = 90. 
#     nstrike = int( length / size_subfault )
#     ndip = int( width / size_subfault )
 
    print("nstrike, ndip: %d %d" % (nstrike,ndip) )
 
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
     
    # Then make the decollement geometry
 
    width_decollement = 100.
    dip_decollement = 0.1
 
    # load the geometry as a np array
 
    np_geometry = np.load('geometry/'+expt+'_geometry.npy')
    print('-- n_dis before adding decollement ', np_geometry.shape[0] )
 
 
    # the top left corner should be with index 0
    # we append it to the end of the file to merge the two geometry
     
    [rdis_long,rdis_lat,rdis_depth,rdis_length,rdis_width,rdis_area,ratio_rdis_tdis,strike,dip,\
     centroid_long,centroid_lat,centroid_depth,\
     tdis_long1,tdis_lat1,tdis_depth1,tdis_long2,tdis_lat2,tdis_depth2,tdis_long3,tdis_lat3,tdis_depth3,\
     tdis_area] = np.copy( np_geometry[0,:] )
 
    # we correctly fill the decollement geometry parameters
     
    # first import the small subfault as a dislocation
    dislocation_lp = Dislocation(0, rdis_long, rdis_lat, rdis_depth, strike, dip, rdis_length, rdis_width, rdis_area, 0., 0.)
 
    # get the corners, the one of interest here is c4
    [c1,c2,c3,c4] = dislocation_lp.corners(coor_type='geo')
    print(c1,c2,c3,c4)
    # for the decollement, we have
    [ rdis_long, rdis_lat, rdis_depth ] = c4
    rdis_length = length
    rdis_width = width_decollement
    # strike decollement = strike subfault so we do not do anything
    dip = dip_decollement
    rdis_area = rdis_length * rdis_width
    # creates a new subfault object and get the centroid
 
#    dislocation_lp_decollement = Dislocation(0, rdis_long, rdis_lat, rdis_depth, strike, dip, length, width, rdis_area, 0., 0.)
    dislocation_lp_decollement = Dislocation(0, -78.50870  ,  0.26632     ,   8.55, strike, dip, length, width, rdis_area, 0., 0.)
     
    centroid = dislocation_lp_decollement.centroid(coor_type='geo')
    [centroid_long,centroid_lat,centroid_depth] = centroid
     
    # append
    decollement_subfault_line = [rdis_long,rdis_lat,rdis_depth,rdis_length,rdis_width,rdis_area,ratio_rdis_tdis,strike,dip,\
     centroid_long,centroid_lat,centroid_depth,\
     tdis_long1,tdis_lat1,tdis_depth1,tdis_long2,tdis_lat2,tdis_depth2,tdis_long3,tdis_lat3,tdis_depth3,\
     tdis_area]
 
    # shallow Ilalo fault
 
#     lon0 = -78.27
#     lat0 = -0.1
#     depth = .0 / 111.
#     strike = 195.
#     dip = 10.
#     length = 40.
#     width  = 20.
#     size_subfault = 2. 
# 
#     nstrike = int( length / size_subfault )
#     ndip = int( width / size_subfault )
# 
#     cmd = 'cd geometry ; pyeq_make_rectangular_fault.py '\
#           + ' -corner /'+str(lon0)+'/'+str(lat0)+'/'+str(depth)\
#           + ' -length '+str(length)\
#           + ' -width '+ str(width)\
#           + ' -strike '+ str(strike)\
#           + ' -dip '+str(dip)\
#           + ' -e ' + 'shallow'\
#           + ' --nstrike '+ str(nstrike)\
#           + ' --ndip '+str(ndip)\
#           + ' --verbose'
#     
#     print('-- running: ',cmd)
#     syscmd.getstatusoutput( cmd , verbose=False )
#
#    # merge the shallow geometry
#    
#    np_geometry_shallow = np.load('geometry/'+'shallow'+'_geometry.npy')
#    np_geometry = np.append(np_geometry,np_geometry_shallow,axis=0)
 
    # let's append the decollement at the end     
    np_geometry = np.append(np_geometry,np.array( decollement_subfault_line ).reshape(1,-1),axis=0)
 
    print('-- n_dis after  adding decollement and shallow ', np_geometry.shape[0] )
     
#     # we correctly fill the length, width and area
# 
#     np_geometry[-1,3] =  length
#     np_geometry[-1,4] =  width_decollement
#     np_geometry[-1,5] =  length * width_decollement
#     np_geometry[-1,8] =  dip_decollement  
# 
#     # now we change the centroid
#     
#     np_geometry[-1,9]  = -78.85
#     np_geometry[-1,10] = -0.1
#     np_geometry[-1,11] = np_geometry[-1,2] # depth
 
    print("-- decollement depth: ",  rdis_depth )
 
    # we rewrite the geometry
     
     
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
# DECOLLEMENT + 3 SEGMENTS GEOMETRY
###############################################################################

    # Decollement
    
    lon0 = -78.43
    lat0 = 0.6
    depth = 13.
    strike = 200.
    dip = .1
    length = 200.
    width  =  200.
    size_subfault = 1000 

    nstrike = int( length / size_subfault )
    ndip = int( width / size_subfault )

    print("- Making decollement" )
    print("nstrike, ndip: %d %d" % (nstrike,ndip) )

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

    np_decollement = np.load('geometry/'+expt+'_geometry.npy')
    
    # ramp central segment Quito
    
    lon0 = -78.4
    lat0 =  -0.10
    depth = 0.
    strike = 200.
    dip = 25.
    length = 38.
    width  =  33.
    size_subfault = 2 

    nstrike = int( length / size_subfault )
    ndip = int( width / size_subfault )

    print("- Making ramp central Quito" )
    print("nstrike, ndip: %d %d" % (nstrike,ndip) )

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

    np_ramp_central = np.load('geometry/'+expt+'_geometry.npy')

    # remove some elements
    long_min_central = -78.42
    lat_max_central = -0.10

#    lindex = np.where( (np_ramp_central[:,0] < long_min_central) &  (np_ramp_central[:,1] < lat_max_central) )
    lindex = np.where( (np_ramp_central[:,0] < long_min_central) )

    np_ramp_central = np_ramp_central[lindex]

    depth_ramp_central = np.max(np.fabs(np_ramp_central[:,2]))

    print('-- geometry central ', np_ramp_central.shape)
    
    # ramp north segment Quito
    
    lon0 = -78.35
    lat0 =   0.25
    depth = 0.
    strike = 200.
    dip = 25.
    length = 39.
    width  =  23.
    size_subfault = 2 

    nstrike = int( length / size_subfault )
    ndip = int( width / size_subfault )

    print("- Making ramp north Quito" )
    print("nstrike, ndip: %d %d" % (nstrike,ndip) )

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

    np_ramp_north = np.load('geometry/'+expt+'_geometry.npy')

    long_min_north = -78.43
    lindex = np.where( (np_ramp_north[:,0] < long_min_north ) )

    np_ramp_north = np_ramp_north[lindex]

    print('-- geometry north ', np_ramp_north.shape)
    
    # ramp south segment Quito
    
    lon0 = -78.58
    lat0 =   -0.4
    depth = 0.
    strike = 200.
    dip = 25.
    length = 35.
    width  =  25.
    size_subfault = 2 

    nstrike = int( length / size_subfault )
    ndip = int( width / size_subfault )

    print("- Making ramp south Quito" )
    print("nstrike, ndip: %d %d" % (nstrike,ndip) )

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

    np_ramp_south = np.load('geometry/'+expt+'_geometry.npy')

    print('-- geometry south ', np_ramp_south.shape)
    
    # merge geometry
    
    print('-- merging geometry')
    
    np_geometry = np.vstack( (np_ramp_north , np_ramp_central) )
    np_geometry = np.vstack( (np_geometry , np_ramp_south) )
    
#    np_geometry = np.append(np_geometry,np_decollement,axis=0)
    np_geometry = np.insert(np_geometry,0,np_decollement,axis=0)
    
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
    
    
###################################################################################
    # Decollement
    
    lon0 = -78.43
    lat0 = 0.6
    depth = 13.
    strike = 200.
    dip = .1
    length = 200.
    width  =  200.
    size_subfault = 1000 

    nstrike = int( length / size_subfault )
    ndip = int( width / size_subfault )

    print("- Making decollement at depth: %1.lf " % (depth) )
    
    print("nstrike, ndip: %d %d" % (nstrike,ndip) )

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

    np_decollement = np.load('geometry/'+expt+'_geometry.npy')

