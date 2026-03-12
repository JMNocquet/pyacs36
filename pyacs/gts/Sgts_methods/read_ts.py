  ###################################################################
def read_ts(self,ts_dir='.', verbose=True, name_filter='', add_key='', sites=[],lexclude=[], type = None, xyz=True):
###################################################################
    """Read time series from directory; format is auto-detected or set by type.

    Supported formats: pos, kenv, mb_file, cats, txyz (pyacs), track (NEU), pride, pck.

    Parameters
    ----------
    ts_dir : str, optional
        Directory of time series files. Default is '.'.
    verbose : bool, optional
        If True, print progress. Default is True.
    name_filter : str, optional
        Filter on time series name ('*name_filter*').
    add_key : str, optional
        String to add before site code.
    sites : list, optional
        If non-empty, only these site codes are read.
    lexclude : list, optional
        Site codes to exclude from reading.
    type : str, optional
        Format: 'pos', 'kenv', 'mb_file', 'cats', 'txyz', 'track', 'pride', 'pck'.
    xyz : bool, optional
        For pos files, read XYZ rather than dNEU. Default is True.

    Returns
    -------
    None
        Time series are appended to this Sgts instance (in-place).
    """
    # import
    import numpy as np
    import sys
    from pyacs.gts.Gts import Gts
    import os

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    #######################################################################
    # CASE PCK FILE HAS BEEN PROVIDED
    # DO NOT ATTEMPT TO READ OTHER FILES
    if os.path.splitext(ts_dir)[1] == '.pck':
        if os.path.isfile(ts_dir):
            lGts={}
            # read Sgts from pickle
            import pickle

            with open(ts_dir, "rb") as f:
                ts = pickle.load(f)
            f.close()

            # Normalize the code to string
            MESSAGE("Normalizing GPS site codes in pickle file")
            try:
                ts.__dict__ = {str(k).strip(): v for k, v in ts.__dict__.items()}
            except Exception as e:
                ERROR(f"Failed to normalize site codes: {str(e)}", exit=True)

            for code in ts.lcode():
                lGts[code] = ts.__dict__[code]
            MESSAGE("Sgts read from PYACS pck file: %s" % (ts_dir))

            for gts in list(lGts.values()):
                try:
                    gts.code = str(gts.code).strip()
                    if len(gts.code) < 4:
                        WARNING(f"Site code {gts.code} is shorter than 4 characters")
                except Exception as e:
                    ERROR(f"Failed to normalize site code for Gts object: {str(e)}", exit=True)
                    
                if gts.code not in lexclude:
                    VERBOSE("adding %s " % gts.code)
                    self.append(gts)

            VERBOSE("read %d time series from %s " % (len(self.lcode()), ts_dir))
            return
        else:
            if not os.path.exists(ts_dir):
                ERROR("File %s does not exist." % ts_dir, exit=True)
            else:
                ERROR("File %s is not a pickle file" % ts_dir, exit=True)
    VERBOSE("Reading time series in directory %s " % ts_dir)
    
    #######################################################################
    def __pride_pos(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
         
        import glob
        import os
         
        lGts={}
        lfile = glob.glob(ts_dir+'/'+'pos_*'+name_filter)
        lcode = []
        for ifile in lfile:
            if not os.path.isfile( ifile ): continue
            code = ifile[-4:].upper()
            if code not in lcode:
                lcode.append( code )
             
        for code in lcode:
            VERBOSE("Reading pride pos files for: %s " % code )
            lGts[code+add_key]=Gts(code=code)
            lGts[code+add_key].read_pride_pos(tsdir=ts_dir , verbose=verbose )
    
        return(lGts)
    
    #######################################################################
    def __pride(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
        import glob
        import os
         
        lGts={}
        lfile = glob.glob(ts_dir+'/'+'kin_*'+name_filter)
        lcode = []
        for ifile in lfile:
            if not os.path.isfile( ifile ): continue
            code = ifile[-4:].upper()
            if code not in lcode:
                lcode.append( code )
             
        for code in lcode:
            VERBOSE("Reading pride pos files for: %s " % code )
            lGts[code+add_key]=Gts(code=code)
            lGts[code+add_key].read_pride(tsdir=ts_dir , verbose=verbose )
    
        return(lGts)
    
    #######################################################################
    def __tdp(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
        import glob
        import os
         
        lGts={}
        ldat=sorted(glob.glob(ts_dir+'/'+'*'+name_filter+'*.tdp'))
        for fdat in ldat:
            if not os.path.isfile( fdat ): continue
            mb_file=os.path.basename(fdat)[0:4]
            code=mb_file[0:4]
            if code not in sites and sites !=[]:continue
            VERBOSE("Reading %s " % mb_file )
            lGts[code+add_key]=Gts(code=code)
            lGts[code+add_key].read_tdp(ifile=fdat,gmt=False)
    
        return(lGts)
    
    #######################################################################
    def __sol(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
        import glob
        import os
         
        lGts={}
        lsol=sorted(glob.glob(ts_dir+'/'+'*'+name_filter+'*.sol'))
        for fsol in lsol:
            if not os.path.isfile( fsol ): continue
            mb_file=os.path.basename(fsol)[0:4]
            code=mb_file[0:4]
            if code not in sites and sites !=[]:continue
            VERBOSE("Reading %s " % mb_file )
            lGts[code+add_key]=Gts(code=code)
            lGts[code+add_key].read_sol(ifile=fsol)
    
        return(lGts)


    #######################################################################
    def __mb_files(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
        import glob
        import os
        lGts={}
        ldat=sorted(glob.glob(ts_dir+'/'+'*'+name_filter+'*.dat1'))
        for fdat in ldat:
            if not os.path.isfile( fdat ): continue
            mb_file=fdat[0:-1]
            code=mb_file[-12:-8]
            if code not in sites and sites !=[]:continue
            VERBOSE("Reading %s " % mb_file)
            lGts[code+add_key]=Gts(code=code)
            lGts[code+add_key].read_mb_file(ifile=mb_file,gmt=False)
    
        return(lGts)
     
    
    #######################################################################
    def __kenv(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
        import glob
        import os
         
        lGts={}
        lkenv=sorted(glob.glob(ts_dir+'/*.kenv'))
        for fkenv in lkenv:
            if not os.path.isfile( fkenv ): continue
            code=fkenv[-20:-9]
            VERBOSE("Reading %s " % fkenv )
            lGts[code+add_key]=Gts(code=code)
            try:
                lGts[code+add_key].read_kenv(fkenv,date_type='jd')
            except:
                ERROR("Could not read %s" % fkenv)
                del lGts[code+add_key]
                continue
        return(lGts)
    
    #######################################################################
    def __cats(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
        import glob
        import os
        lGts={}
        lcats=sorted(glob.glob(ts_dir+'/cats*.dat'))
        for cats in lcats:
            if not os.path.isfile( cats ): continue
            code=cats.split('cats_')[-1][:4].upper()
            VERBOSE("Reading %s" % cats)
            lGts[code+add_key]=Gts(code=code)
            try:
                lGts[code+add_key].read_cats_file(ifile=cats,gmt=False)
            except:
                ERROR("Could not read %s" % cats )
                del lGts[code+add_key]
                continue
        return(lGts)
    
    #######################################################################
    def __pos(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose,xyz=True):
        import glob
        import os


        lGts={}
        # change if a single pos file is provided
        if os.path.isfile(ts_dir) and os.path.splitext(ts_dir)[1] == '.pos':
            lpos = [ ts_dir ]
        else:
            lpos=sorted(glob.glob(ts_dir+'/????*.pos'))

        for pos_w_path in lpos:
            if not os.path.isfile( pos_w_path): continue
            _drive, path_and_file = os.path.splitdrive(pos_w_path)
            _path, pos = os.path.split(path_and_file)
            code=pos[:4]
            if code not in sites and sites !=[]:continue
            VERBOSE("Reading %s" % pos)
            wts = Gts().read_pos(tsfile=pos_w_path,xyz=xyz)
            lGts[wts.code+add_key]=wts
            #lGts[code+add_key]=Gts(code=code)
            #lGts[code+add_key].read_pos(tsfile=pos_w_path,xyz=xyz)
        return(lGts)
    
    #######################################################################
    def __ts_xyz(ts_dir='.',name_filter=None, verbose=False):
         
        import glob
        import os
         
        lGts={}
        l_txyz=glob.glob(ts_dir+'/*txyz.dat')
         
             
        ##############################################
        def __read_txyz(file_name,lGts,verbose=False):
    
            """
            Reads a single txyz.dat
            """
             
            NP_STR_TXYZ=np.genfromtxt(file_name,dtype=str)
            NP_CODE=NP_STR_TXYZ[:,2]
            NP_DATE=np.array(NP_STR_TXYZ[:,1],dtype=float)
            NP_XYZ=np.array(NP_STR_TXYZ[:,4:],dtype=float)
         
            for i in np.arange(NP_CODE.shape[0]):
                code = NP_CODE[i]
                 
                if code in list(lGts.keys()):
                    lGts[code].add_obs_xyz(NP_DATE[i],NP_XYZ[i],in_place=True,check=False)
                else:
                    new_gts=Gts(code=code)
                    new_gts.add_obs_xyz(NP_DATE[i],NP_XYZ[i],in_place=True,check=False)
                    lGts[code]=new_gts
         
            return(lGts)
         
        i=1
        for txyz in l_txyz:
            if not os.path.isfile( txyz ): continue
            VERBOSE("Reading %s #%05d / %05d" % (txyz,i,len(l_txyz)))
            i = i + 1
                 
            lGts=__read_txyz(txyz, lGts,verbose=verbose)
         
        for code in sorted( lGts.keys() ):
             
            gts = lGts[code]
             
            VERBOSE("Generating NEU from XYZ for %s" % gts.code)
            gts.xyz2neu(corr=False)
             
            if not gts.cdata():
                sys.exit()
             
        return(lGts)
         
    
    #######################################################################
    def __track_NEU(ts_dir='.', name_filter=name_filter, add_key=add_key, verbose=verbose):
        import glob
        lGts={}
        ltrack=sorted(glob.glob(ts_dir+'/*.NEU.*'))
        import os
        for track_w_path in ltrack:
            if not os.path.isfile( track_w_path ): continue
            _drive, path_and_file = os.path.splitdrive(track_w_path)
            _path, track = os.path.split(path_and_file)
            code= track.split('.')[-2].upper()
            if code not in sites and sites !=[]:continue
            VERBOSE("Reading %s" % track)
            lGts[code+add_key]=Gts(code=code)
            try:
                lGts[code+add_key].read_track_NEU(tsfile=track_w_path)
            except:
                ERROR("Could not read %s" % track)
                del lGts[code+add_key]
                continue
        return(lGts)

    # END OF SUBROUTINES  #######################################################################

    #######################################################################
    # START OF MAIN
    #######################################################################

    lGts={}
    #######################################################################
    # PCK FORMAT
    #######################################################################
    # Added JMN 22/08/2019
    # Sgts saved as a pickle dump
    if type is None or type=='pck':
        import glob

        if os.path.isfile( self.dir ) and self.dir[-3:]=='pck':
            lpck = [self.dir]
        else:
            lpck = glob.glob( self.dir+'/*.pck' )

        if lpck != []:
            pck = lpck[-1]
            # read Sgts from pickle
            import pickle
            with open( pck, "rb") as f:
                ts = pickle.load(f)
            f.close()
             
            for code in ts.lcode():
                lGts[code] = ts.__dict__[code]
            VERBOSE("Sgts read from PYACS pck file: %s" % (pck) )
        else:
        # no pck file
            VERBOSE("No PYACS pck file found.")
    
    
    #######################################################################
    # PRIDE POS
    #######################################################################
    if type is None or type=='pride_pos':
        lGts_pride_pos = __pride_pos(ts_dir=self.dir,verbose=verbose)
        if lGts_pride_pos=={}:
            VERBOSE("No pride pos files found")
        else:
            lGts=lGts_pride_pos
    
    #######################################################################
    # PRIDE
    #######################################################################
    if type is None or type=='pride':
        lGts_pride = __pride(ts_dir=self.dir,verbose=verbose)
        if lGts_pride=={}:
            VERBOSE("No pride_files found")
        else:
            lGts=lGts_pride
    
    #######################################################################
    # MB_FILE
    #######################################################################
    if type is None or type=='mb_file':
        lGts_mb_files=__mb_files(ts_dir=self.dir,verbose=verbose)
        if lGts_mb_files=={}:
            VERBOSE("No mb_files found")
        else:
            lGts=lGts_mb_files
             
            try:
                self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
            except:
                WARNING("Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for mb_files.")
                WARNING("Some PYACS function will produce errors")
     
    #######################################################################
    # TDP
    #######################################################################
    if type is None or type=='tdp':
        lGts_mb_files=__tdp(ts_dir=self.dir,verbose=verbose)
        if lGts_mb_files=={}:
            VERBOSE("No tdp_files found")
        else:
            lGts=lGts_mb_files
            try:
                self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
            except:
                WARNING("Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for tdp file.")
                WARNING("Some PYACS function will produce errors")

    #######################################################################
    # KENV
    #######################################################################
    if type is None or type=='kenv':
        lGts_kenv=__kenv(ts_dir=self.dir,verbose=verbose)
        if lGts_kenv=={}:
            VERBOSE("No kenv file found")
        else:lGts=lGts_kenv

    #######################################################################
    # SOL
    #######################################################################
    if type is None or type=='sol':
        lGts_sol=__sol(ts_dir=self.dir,verbose=verbose)
        if lGts_sol=={}:
            VERBOSE("No sol file found")
        else:lGts=lGts_sol


    #######################################################################
    # CATS
    #######################################################################
    if type is None or type=='cats':
        lGts_cats=__cats(ts_dir=self.dir,verbose=verbose)
        if lGts_cats=={}:
            VERBOSE("No cats file found")
        else:
            lGts=lGts_cats
            try:
                self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
            except:
                WARNING("Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for cats files")
                WARNING("Some PYACS function will produce errors")

     
    #######################################################################
    # POS FORMAT
    #######################################################################
    if type is None or type=='pos':
        lGts_pos=__pos(ts_dir=self.dir,verbose=verbose,xyz=xyz)
        if lGts_pos=={}:
            VERBOSE("No Gamit/Globk pos file found")
        else:lGts=lGts_pos
    
    #######################################################################
    # TXYZ FORMAT
    #######################################################################
    if type is None or type=='txyz':
        lGts_txyz=__ts_xyz(ts_dir=self.dir,verbose=verbose)
        if lGts_txyz=={}:
            VERBOSE("No pyacs t_xyz file found")
        else:lGts=lGts_txyz
    
    #######################################################################
    # TRACK FORMAT
    #######################################################################
    if type is None or type=='track':
        lGts_track=__track_NEU(ts_dir=self.dir,verbose=verbose)
        if lGts_track=={}:
            VERBOSE("No Gamit/Globk track NEU file found")
        else:lGts=lGts_track
    
    
    for gts in list(lGts.values()):
        if gts.code not in lexclude:
            VERBOSE("adding %s " % gts.code)
            self.append(gts)
    
    VERBOSE("read %d time series from %s " % (len(self.lcode()), ts_dir))
    
    
    
