"""
Reads old GAMIT/GLOBK time series mb_files
"""

###################################################################
## Loads Gts from GAMIT/GLOBK mb_file format for time series
###################################################################
def read_mb_file(self,idir='.',ifile=None, gmt=True, verbose=False):
    """
    Read GAMIT/GLOBK mb_files in a directory and actually loads the time series
    
    """
    
    import numpy as np
    import os
    
    if gmt==True:
        gmt_file=idir+'/../maps/en_velo.gmt'
    if  isinstance(gmt,str):
        gmt_file=gmt
    
    if gmt != False:
        self.read_lon_lat(gmt_file,verbose=verbose)
    
    if ifile is None:
        mb_file_basename= idir + '/mb_'+self.code+'_GPS.dat'
    else:
        mb_file_basename=ifile
    
    data_NEU = []
    for i in range(1,4):
        mb_file = mb_file_basename + str(i)

        # file
        self.ifile=os.path.abspath(mb_file)
        
        data=np.genfromtxt(mb_file,skip_header=4)
        
    # reshape to ensure a 2D array
        if len(data.shape)==1:
            data=data.reshape((1,data.shape[0]))
        


        data_NEU.append(data)

    if data_NEU[0].shape == data_NEU[1].shape == data_NEU[2].shape:
        self.data=np.zeros((data_NEU[0].shape[0],7))
        self.data[:,0]=data_NEU[0][:,0]
        self.data[:,1]=data_NEU[0][:,1]#*to_mm
        self.data[:,2]=data_NEU[1][:,1]#*to_mm
        self.data[:,3]=data_NEU[2][:,1]#*to_mm

        self.data[:,4]=data_NEU[0][:,2]#*to_mm
        self.data[:,5]=data_NEU[1][:,2]#*to_mm
        self.data[:,6]=data_NEU[2][:,2]#*to_mm

    else: 
        print("!!! Error reading ",mb_file_basename," :*dat1, *dat2, *dat3 do not have the same length")
        self.data = None

###################################################################
def write_mb_file(self,idir,add_key=''):
###################################################################
    dat1=open((idir + '/mb_' + self.code.upper() + add_key+'_GPS.dat1'),'w+')
    dat1.write("Pyacs Analysis\n" + self.code.upper() + add_key+"_GPS to N Solution  1\n\n\n");
    dat2=open((idir + '/mb_' + self.code.upper() + add_key+'_GPS.dat2'),'w+')
    dat2.write("Pyacs Analysis\n" + self.code.upper() + add_key+"_GPS to E Solution  1\n\n\n");
    dat3=open((idir + '/mb_' + self.code.upper() + add_key+'_GPS.dat3'),'w+')
    dat3.write("Pyacs Analysis\n" + self.code.upper() + add_key+"_GPS to U Solution  1\n\n\n");


    ndate=self.data.shape[0] 
    for i in range(ndate):
        ### for North composant
        dat1.write(' %10.5lf   ' % self.data[i,0]), dat1.write('%10.5lf' % (self.data[i,1])), dat1.write(' %10.5lf\n' % (self.data[i,4]))
        ### for East composant
        dat2.write(' %10.5lf   ' % self.data[i,0]), dat2.write('%10.5lf' % (self.data[i,2])), dat2.write(' %10.5lf\n' % (self.data[i,5]))
        ### for Up composant
        dat3.write(' %10.5lf   ' % self.data[i,0]), dat3.write('%10.5lf' % (self.data[i,3])), dat3.write(' %10.5lf\n' % (self.data[i,6]))
    
    dat1.close
    dat2.close
    dat3.close
    return(self)
