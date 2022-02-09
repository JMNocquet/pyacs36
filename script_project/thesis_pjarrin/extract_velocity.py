#!/usr/bin/env python

# pos files to be analyzed
pos_dir = '/Users/nocquet/projets/2018/soam_proc/run_pjarrin_12_12_2018/pos'

# outdir
outdir = '/Users/nocquet/projets/2018/soam_proc/pyacs_analysis_velocity_field'

# outfiles
cgps_file ='vel_cgps.dat'
sgps_file ='vel_sgps.dat'
ssgps_file ='vel_ssgps.dat'

from pyacs.gts.Sgts import Sgts

ts = Sgts(ts_dir=pos_dir)

# select sites impacted by the 2016/04/16 EQ

bounds_eq_2016 = [-83,-74,-3.5, 2.5]
period = [1990., 2016.29]

ts_eq_2016 = ts.sel_rectangle(bounds_eq_2016, verbose=True)

# loop on sites

for gts in ts.lGts():
    print("-- Processing " , gts.code)

    # remove post EQ data if it is in the selected area

    if gts.code in ts_eq_2016.lcode():
        print("-- EQ 2016")
        gts = gts.extract_periods(period) 

    if gts.data is None:
        continue

    # detrend_median pour cgps data
    detrended = gts.detrend_median()

    if detrended is not None:
        print("-- saving CGPS vel to " , cgps_file )
        detrended.save_velocity(cgps_file)
        continue
    
    # detrend for sgps data
    detrended = gts.detrend_median(delta_day=0)
    if detrended is not None:
        print("-- saving SGPS vel to " , sgps_file )
        detrended.save_velocity(sgps_file)
        continue
    
    # detrend for ssgps data (only two measurements)
    detrended = gts.detrend()
    if detrended is not None:
        print("-- saving SGPS vel to " , ssgps_file )
        detrended.save_velocity(ssgps_file)
        continue
    