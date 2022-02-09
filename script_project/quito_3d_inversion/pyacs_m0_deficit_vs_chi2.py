#!/usr/bin/env python

mdir = '/Users/nocquet/projets/2018/quito/models_03_12_2018/model'

from glob import glob
import numpy as np
import re

mu=3.E10

lmodel_vh_45 = glob(mdir+'/*vh_04*sigma_01_dc_4*')

lm0=[]
lchi2=[]
lcoupling=[]
for model_dir in lmodel_vh_45:
    print("-- analyzing %s" % model_dir)
    sol_slip = glob(model_dir+'/*sol_slip.dat')[-1]
    model = np.genfromtxt(sol_slip)
    # remove the decollement
    model=np.delete(model, 0, axis=0)
    # reference velocity
    dip = model[0,8]
    print("--dip ", dip)
    v = 4.5 / np.cos(np.radians(dip))
    print("-- v ", v)
    
    area = np.sum( model[:,5] )
    
    max_moment = area * 1.E6 * 4.5 * 1.E-3 * mu
    actual_moment = area * 1.E6 * np.mean( model[:,-1] )  * 1.E-3 * mu

    coupling = (max_moment - actual_moment)/max_moment * 100.
    
    print("-- area ", area)
    print('-- maximum moment :', max_moment)
    print('-- actual  moment :', actual_moment)
    print('-- moment  deficit:', max_moment - actual_moment)
    print("-- average coupling %.1lf " % coupling )
    
    sum = glob(model_dir+'/*sum.dat')[-1]
    
    with open(sum) as f:
        content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        content = [x.strip() for x in content] 
    chi2  = list(filter(lambda x: re.search(r'chi2 all obs', x), content))[0].split(':')[-1].lstrip()
    print("-- chi2" , chi2 )
    lm0.append(max_moment - actual_moment)
    lchi2.append(chi2)
    lcoupling.append(coupling)

d=np.array( [lm0,lchi2,lcoupling],dtype=float).T
np.savetxt('moment_deficit_vs_chi2_vh_45.dat', d, fmt="%5.3E %5.3E %4.1lf")