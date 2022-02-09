from pyeq.lib.okada.okada import okada

import numpy as np

llambda = 3.2074E+10  
mu = 3.9701E+10

SOURCE = False

if SOURCE:
    
    # source
    
    slip = 1.
    area = 1. * 1 
    
    dislocations = np.array ( slip * area ).reshape(-1,1)
    
    dips = np.array ( 45. ).reshape(-1,1) 
    rakes = np.array( 90. ).reshape(-1,1)
    
    xs = np.array( 0. ).reshape(-1,1)
    ys = np.array( 0. ).reshape(-1,1)
    zs = np.array( 1. ).reshape(-1,1)
    
    lengths = np.array( 0. ).reshape(-1,1)
    widths = np.array(  0. ).reshape(-1,1)
    strikes = np.array( 0. ).reshape(-1,1)
    
    # gps
    
    xrec  = np.array( 0. ).reshape(-1,1)
    yrec  = np.array( 1. ).reshape(-1,1)
    zrec0 = 0.
    
    # run okada
    (disp,tilt,strain)=okada(llambda,mu,dislocations, xs,ys,zs,lengths,widths,strikes,dips,rakes, xrec,yrec,zrec0)
#    print("disp:   %7.4E %7.4E %7.4E" % tuple(disp.squeeze()) )
#    print("tilt:   %7.4E %7.4E " % tuple(tilt.squeeze()) )
#    print("strain: %7.4E %7.4E %7.4E %7.4E %7.4E %7.4E" % tuple(strain.squeeze()) )


# Test dc3d

xrec = np.array([-8.E3 , -20.E3, -32.E3, -150.E3 ])
yrec = np.array([100.E3 , 100.E3, 100.E3, 100.E3 ])

#xrec = np.array([[-8.E3  ]])
#yrec = np.array([[100.E3 ]])

zrec0 = 0.

slips = np.array ( 5. ).reshape(-1,1)
xs = np.array( 0. ).reshape(-1,1)
ys = np.array( 0. ).reshape(-1,1)
zs = np.array( 30.E3 ).reshape(-1,1)

lengths = np.array( 200.E3 ).reshape(-1,1)
widths = np.array(  50.E3 ).reshape(-1,1)
strikes = np.array( 90. ).reshape(-1,1)
dips = np.array ( 12. ).reshape(-1,1) 
rakes = np.array( 90. ).reshape(-1,1)

(disp,tilt,strain)=okada(llambda,mu,slips, xs,ys,zs,lengths,widths,strikes,dips,rakes, xrec,yrec,zrec0)

#print('disp')
#print(disp)

#print('tilt')
#print(tilt)

#print('strain')
#print(strain[:,4:7])

#print("disp:   %7.4E %7.4E %7.4E" % tuple(disp.squeeze()) )
#print("tilt:   %7.4E %7.4E " % tuple(tilt.squeeze()) )
#print("strain: %7.4E %7.4E %7.4E %7.4E %7.4E %7.4E" % tuple(strain.squeeze()) )

