
"""
Pyacs noise module
This is a wrapper to Globk tsfit and the CATS Software
Reference
Herring, T. A., King, R. W., Floyd, M. A. & McClusky, S. C. (2015). GAMIT/GLOBK Reference manual, 10.6. 
Williams, S. D. (2008). CATS: GPS coordinate time series analysis software. GPS solutions, 12(2), 147-153.
"""

import numpy as np
import os, sys

###############################################################################
def wrms(self):
###############################################################################
    """
    Return the wrms
    :return wrms: return(np.array([wrms_n,wrms_e,wrms_up]))
    """

    data_pos   = self.data[:,1:4]
    data_sigma = self.data[:,4:7]
    
    wrms = np.sqrt( np.sum( (data_pos / data_sigma)**2 , axis=0 ) / np.sum( 1. / data_sigma**2 , axis=0 ) )
    
    return(wrms)

    
###############################################################################
def realistic_sigma(self,option='tsfit',in_place=False,verbose=False):
###############################################################################
    """
    Calculates realistic sigmas on velocity components
    :param option:
        - tsfit: globk T. Herring realistic sigma
        - cats_pl: CATS estimates with noise type estimated (i.e. --model=pl:)
        - cats_seasonal_pl: CATS estimates with seasonal terms and noise type estimated (i.e. --model=pl: --sinusoid=1y1)
        - cats_flicker: CATS estimates assuming flicker noise (i.e. --model=pl:k-1)
        - cats_seasonal_flicker: CATS estimates with seasonal terms and assuming flicker noise (i.e. --model=pl:k-1 --sinusoid=1y1)
    """

    if option=='cats_pl':
        return(self.sigma_cats(in_place=in_place,verbose=verbose,k='',seasonal=''))
    if option=='cats_seasonal_pl':
        return(self.sigma_cats(in_place=in_place,verbose=verbose,k='',seasonal='True'))
    if option=='cats_flicker':
        return(self.sigma_cats(in_place=in_place,verbose=verbose,k='k-1',seasonal=''))
    if option=='cats_seasonal_flicker':
        return(self.sigma_cats(in_place=in_place,verbose=verbose,k='k-1',seasonal='True'))

    if option=='tsfit':
        return(self.sigma_vel_tsfit(in_place=in_place,verbose=verbose))


###############################################################################
def sigma_cats(self,in_place=False,verbose=False,k='k-1',seasonal=''):
###############################################################################
    """
    runs CATS for getting realistic sigma
    """
    # setup_save working directory    
    import shutil

    name_dir=('tmp_cats_%s' % self.code)
    
    try:
        shutil.rmtree(name_dir)
    except:
        pass
    
    try:
        os.mkdir(name_dir)
    except:
        print('!!! ERROR could not create directory ',name_dir)
        sys.exit()


    if verbose:
        print("-- dir ",name_dir,' created')
    
    os.chdir(name_dir)

    # write cats file

    self.write_cats('.')

    # runs cats
    
    import subprocess
    
    cats_file=("cats_%s.dat" % self.code)
    
    if seasonal !='':
        seasonal_opt='--sinusoid=1y1'
    else:
        seasonal_opt=''
    
    cmd=('cats %s --model=pl:%s %s -o cats_output.dat' % (cats_file,k,seasonal_opt)) 
    
    if verbose:
        print('  -- Running ',cmd)
        cmd_output=subprocess.getstatusoutput(cmd)
    else:
        subprocess.getstatusoutput(cmd)

    # read results




  
    try:
        f=open('cats_output.dat')
        f.close()
    except:
        print('!!! ERROR: problem running CATS')
        return('ERROR')
    
    
    if 'k' not in k:
        (ke,kn,ku)=get_spectral_index('cats_output.dat')
        print(("-- Noise spectral index : N=%.1lf E=%.1lf U=%.1lf" % (ke,kn,ku)))
       
    (Ve,Vn,Vu,SVe,SVn,SVu)=get_vel_and_sigma('cats_output.dat')
    
    os.chdir('..')
    
    # verbose
    
    if verbose:
        print((" --  Ve %8.2lf  Vn %8.2lf  Vu %8.2lf" % (Ve,Vn,Vu)))
        print((" -- SVe %8.2lf SVn %8.2lf SVu %8.2lf" % (SVe,SVn,SVu)))
    
    # populates and return new time series
    
    new_gts=self.copy()

    if new_gts.velocity is None:
        new_gts.velocity=np.zeros(6)
        new_gts.velocity[0]=Vn/1000.0
        new_gts.velocity[1]=Ve/1000.0
        new_gts.velocity[2]=Vu/1000.0



    new_gts.velocity[3]=SVn/1000.0
    new_gts.velocity[4]=SVe/1000.0
    new_gts.velocity[5]=SVu/1000.0

    # returns    
    if in_place:
            self.velocity=new_gts.velocity
    else:
        return(new_gts)



###############################################################################
def sigma_vel_tsfit(self,in_place=False,verbose=False):
###############################################################################
    """
    runs tsfit for getting realistic sigma
    """

    # setup_save working directory    
    import shutil

    name_dir=('tmp_tsfit_%s' % self.code)
    
    try:
        shutil.rmtree(name_dir)
    except:
        pass
    
    try:
        os.mkdir(name_dir)
    except:
        print('!!! ERROR could not create directory ',name_dir)
        sys.exit()


    if verbose:
        print("-- dir ",name_dir,' created')
    
    os.chdir(name_dir)
    
    # write pos file

    self.write_pos('.')
    
    # creates tsfit.cmd file
    
    tsfit_cmd='tsfit.cmd'
    
    f_tsfit_cmd=tsfit_cmd
    f=open(f_tsfit_cmd,'w')
    f.write(" REAL_SIGMA\n")
    f.write(" VELFILE tsfit_vel.dat\n")
    f.close()
    
    # run tsfit 
    
    import subprocess
    
    pos=("%s.pos" % self.code)
    
    cmd=('tsfit %s sum.dat %s' % (tsfit_cmd,pos)) 
    
    if verbose:
        print('  -- Running ',cmd)
        cmd_output=subprocess.getstatusoutput(cmd)
        for line in cmd_output:
            print(line)
    else:
        subprocess.getstatusoutput(cmd)
    
    
    # reads results
    
    new_velocity=np.genfromtxt('tsfit_vel.dat', skip_header=5, usecols=(2,3,6,7,9,11))
    
    [Ve,Vn,SVe,SVn,Vu,SVu]=new_velocity
    
    os.chdir('..')
    
    # verbose
    
    if verbose:
        print((" --  Ve %8.2lf  Vn %8.2lf  Vu %8.2lf" % (Ve,Vn,Vu)))
        print((" -- SVe %8.2lf SVn %8.2lf SVu %8.2lf" % (SVe,SVn,SVu)))
    
    # populates and return new time series
    
    new_gts=self.copy()
    
    if new_gts.velocity is None:
        new_gts.velocity=np.zeros(6)
        new_gts.velocity[0]=new_velocity[1]/1000.0
        new_gts.velocity[1]=new_velocity[0]/1000.0
        new_gts.velocity[2]=new_velocity[2]/1000.0

    
    new_gts.velocity[3]=new_velocity[2]/1000.0
    new_gts.velocity[4]=new_velocity[3]/1000.0
    new_gts.velocity[5]=new_velocity[5]/1000.0


    # returns    
    if in_place:
            self.velocity=new_gts.velocity
    else:
        return(new_gts)

###############################################################################
def get_spectral_index(cats_file):
###############################################################################
    f=open(cats_file,'r')
    for line in f:
        if ('+NORT' in line) and ('INDEX' in line) and ('+-' in line):
            lline=line.split()
            fin_n=float(lline[3])
        if ('+EAST' in line) and ('INDEX' in line) and ('+-' in line):
            lline=line.split()
            fin_e=float(lline[3])
        if ('+VERT' in line) and ('INDEX' in line) and ('+-' in line):
            lline=line.split()
            fin_u=float(lline[3])
    return(fin_e,fin_n,fin_u)

###############################################################################
def get_vel_and_sigma(cats_file):
###############################################################################
    f=open(cats_file,'r')
    for line in f:
        if ('+NORT' in line) and ('SLOPE' in line) and ('+-' in line):
            lline=line.split()
            v_n=float(lline[3])
            sv_n=float(lline[-1])
        if ('+EAST' in line) and ('SLOPE' in line) and ('+-' in line):
            lline=line.split()
            v_e=float(lline[3])
            sv_e=float(lline[-1])
        if ('+VERT' in line) and ('SLOPE' in line) and ('+-' in line):
            lline=line.split()
            v_u=float(lline[3])
            sv_u=float(lline[-1])
    return(v_e,v_n,v_u,sv_e,sv_n,sv_u)
