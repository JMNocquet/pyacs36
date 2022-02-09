"""
Various routines for pyacs/sol module
"""

###############################################################################
def output_directories(wdir):
###############################################################################
    """
    Creates or reads pyacs output directory
    """

    import os

    path_stat=wdir + '/stat'
    path_res_pos=wdir + '/res_pos'
    path_maps=wdir + '/maps'
    path_pos=wdir + '/pos'
    path_pck=wdir + '/pck'
    path_txyz=wdir + '/tsxy'
    path_tsg=wdir + '/tsg'
    path_prt=wdir + '/glred_prt'
    path_lgred_pos=wdir + '/glred_pos'
    path_log=wdir + '/glred_log'
    path_sinex=wdir + '/glred_sinex'
    path_ssc=wdir + '/glred_ssc'
    
    
    try:
        os.mkdir(wdir)
        os.mkdir(path_stat)
        os.mkdir(path_res_pos)
        os.mkdir(path_maps)
        os.mkdir(path_pos)
        os.mkdir(path_pck)
        os.mkdir(path_txyz)
        os.mkdir(path_tsg)
    
        os.mkdir(path_prt)
        os.mkdir(path_lgred_pos)
        os.mkdir(path_log)
        os.mkdir(path_sinex)
        os.mkdir(path_ssc)
    except:
        pass

    return(path_stat,path_res_pos,path_maps,path_pos,path_pck,path_txyz,path_tsg,path_prt,path_lgred_pos,path_log,path_sinex,path_ssc)



###############################################################################
def write_tsg(lsinex,TSG,ofile):
###############################################################################
    """
    Write Helmert parameters
    """
    
    import numpy as np
    
    f=open(ofile,'a')
    
    f.write('----------------------------------------------------------------------------------------------------------------\n')
    f.write('SINEX                     Epoch        Tx(cm)     Ty(cm)     Tz(cm)      D(ppb)    Rx(mas)    Ry(mas)    Rz(mas)\n')
    f.write('----------------------------------------------------------------------------------------------------------------\n')
    
    
    for i in np.arange(TSG.shape[0]):
        f.write("%12s %11.4lf %11.2f%11.2f%11.2f %11.2f %11.2f%11.2f%11.2f\n" % \
                (lsinex[i].split('/')[-1],TSG[i,0],TSG[i,1]*100.,TSG[i,2]*100.,TSG[i,3]*100.,TSG[i,4]*1000.,TSG[i,5]*1000.,TSG[i,6]*1000.,TSG[i,7]*1000.) )

    f.close()


###############################################################################
def helmert_residuals(RESIDUALS,H_STAT,name,ofile,verbose=True):
###############################################################################
    """
    Print residuals of a Helmert transformation
    """
    import numpy as np

    f=open(ofile,'a')


    f.write('# sinex: %s\n' % name )
    f.write('----------------------------------------------------------------------------\n')
    f.write("n_site: %5d n_obs: %5d n_rejected: %5d %5.1lf %% user_rejection_threshold: %2.1lf\n" % \
            (H_STAT['n_site'],H_STAT['n_obs_raw'],H_STAT['n_obs_rejected'],H_STAT['percent_rejected'],H_STAT['rejection_threshold']))
    f.write('----------------------------------------------------------------------------\n')
    f.write('POINT SOLN  R_East(mm)  R_North(mm)  R_Up(mm)  Outlier_rejected\n')
    f.write('----------------------------------------------------------------------------\n')

    SORTED_RESIDUALS=RESIDUALS[np.argsort(RESIDUALS[:, 0])]
#    SORTED_RESIDUALS=RESIDUALS
    
    min_wrms_2D=1.E6
    min_wrms_2D_site='XXXX'

    min_wrms_3D=1.E6
    min_wrms_3D_site='XXXX'

    max_wrms_2D=0
    max_wrms_2D_site='XXXX'

    max_wrms_3D=0
    max_wrms_3D_site='XXXX'

    
    for i in np.arange(SORTED_RESIDUALS.shape[0]):
        code,soln=SORTED_RESIDUALS[i,:2]
        Re,Rn,Ru=list(map(float,SORTED_RESIDUALS[i,2:5]))
        Rje,Rjn,Rju=list(map(float,SORTED_RESIDUALS[i,5:]))

        SRje='     '
        SRjn='     '
        SRju='     '

        if Rje==1:SRje=' East'
        if Rjn==1:SRjn='North'
        if Rju==1:SRju='   Up'
        
        f.write("%4s %4s %10.2lf  %10.2lf  %10.2lf  %5s %5s %5s\n" % \
                (code,soln,Re,Rn,Ru,SRje,SRjn,SRju))

        if (SRje=='     ') and (SRjn=='     '):
            wrms_2D=np.sqrt(Re**2+Rn**2)
            
            if wrms_2D < min_wrms_2D:
                min_wrms_2D=wrms_2D
                min_wrms_2D_site=code

            if wrms_2D > max_wrms_2D:
                max_wrms_2D=wrms_2D
                max_wrms_2D_site=code


        if (SRje=='     ') and (SRjn=='     ') and (SRju=='     '):
            wrms_3D=np.sqrt(Re**2+Rn**2+Ru**2)
            
            if wrms_3D < min_wrms_3D:
                min_wrms_3D=wrms_3D
                min_wrms_3D_site=code

            if wrms_3D > max_wrms_3D:
                max_wrms_3D=wrms_3D
                max_wrms_3D_site=code



    f.write('----------------------------------------------------------------------------\n')
    f.write("Postfit L1 norm median: %10.2lf %10.2lf %10.2lf \n" % tuple(H_STAT['median_L1']))
    f.write("Postfit L2 norm wrms  : %10.2lf %10.2lf %10.2lf \n" % tuple(H_STAT['wrms']))
    f.write('----------------------------------------------------------------------------\n')
    f.write("Best  fit 2D: %s %10.2lf     Best  fit 3D: %s  %10.2lf \n" % (min_wrms_2D_site,min_wrms_2D,min_wrms_3D_site,min_wrms_2D))
    f.write("Worst fit 2D: %s %10.2lf     Worst fit 3D: %s  %10.2lf \n" % (max_wrms_2D_site,max_wrms_2D,max_wrms_3D_site,max_wrms_3D))
    f.write('----------------------------------------------------------------------------\n')

    f.close()

