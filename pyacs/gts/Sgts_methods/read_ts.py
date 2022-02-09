###################################################################
def read_ts(self,ts_dir='.', verbose=True, name_filter='', add_key='', sites=[],lexclude=[], type = None, xyz=True):
###################################################################
    """
    Reads time series, trying to guess the format. Current time series format supported are: pos, kenv, mb_file, cats, txyz (pyacs), track (NEU format for high rate)
     
    :param ts_dir: directory of time series files
    :param name_filter: string used to filter time series name '*name_filter*'
    :param add_key: adds a string before site code
    :param sites: list of site codes to be read. Any other will be discarded.
    :param lexclude: list of sites to be excluded from reading
    :param type: specifies the format of the time series files. Choose among ['pos', 'kenv', 'mb_file', 'cats', 'txyz', 'track' , 'pride','pck']
    :param xyz: for pos files, reads the XYZ coordinates rather than dNEU. This is the default.
     
    :return: an Sgts instance.
    """
    # import
    import numpy as np
    import sys
    from pyacs.gts.Gts import Gts
    
    if verbose:
        print('-- Reading directory: ',ts_dir)
    
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
            if verbose:
                print("-- Reading pride pos files for: %s " % code )
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
            if verbose:
                print("-- Reading pride files for: %s " % code )
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
            mb_file=fdat.split('/')[-1][0:4]
            code=mb_file[0:4]
            if code not in sites and sites !=[]:continue
            if verbose:
                print("-- Reading ", mb_file)
            lGts[code+add_key]=Gts(code=code)
            lGts[code+add_key].read_tdp(ifile=fdat,gmt=False)
    
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
            if verbose:
                print("-- Reading ", mb_file)
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
            if verbose:print("-- Reading ", fkenv)
            lGts[code+add_key]=Gts(code=code)
            try:
                lGts[code+add_key].read_kenv(fkenv,date_type='jd')
            except:
                print("!!! Error. Could not read ",fkenv)
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
            if verbose:print("-- Reading ", cats)
            lGts[code+add_key]=Gts(code=code)
            try:
                lGts[code+add_key].read_cats_file(ifile=cats,gmt=False)
            except:
                print("!!! Error. Could not read ",cats)
                del lGts[code+add_key]
                continue
        return(lGts)
    
    #######################################################################
    def __pos(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose,xyz=True):
        import glob
        lGts={}
        lpos=sorted(glob.glob(ts_dir+'/????*.pos'))
        import os
        for pos_w_path in lpos:
            if not os.path.isfile( pos_w_path): continue
            _drive, path_and_file = os.path.splitdrive(pos_w_path)
            _path, pos = os.path.split(path_and_file)
            code=pos[:4]
            if code not in sites and sites !=[]:continue
            if verbose:print("-- Reading ", pos)
            lGts[code+add_key]=Gts(code=code)
            lGts[code+add_key].read_pos(tsfile=pos_w_path,xyz=xyz)
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
            if verbose:
                print('-- Reading ',txyz,i,'/',len(l_txyz))
                i = i + 1
                 
            lGts=__read_txyz(txyz, lGts,verbose=verbose)
         
        for code in sorted( lGts.keys() ):
             
            gts = lGts[code]
             
            if verbose:
                print('-- Generating NEU from XYZ for ',gts.code)
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
            if verbose:print("-- Reading ", track)
            lGts[code+add_key]=Gts(code=code)
            try:
                lGts[code+add_key].read_track_NEU(tsfile=track_w_path)
            except:
                print("!!! Error. Could not read ",track)
                del lGts[code+add_key]
                continue
        return(lGts)
    
    
    # START
    
    lGts={}
    #######################################################################
    # PCK FORMAT
    #######################################################################
    # Added JMN 22/08/2019
    # Sgts saved as a pickle dump
    if type is None or type=='pck':
        import glob
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
            if verbose:
                print("-- Sgts read from PYACS pck file: %s" % (pck) )
        else:
        # no pck file    
            if verbose:print("-- No PYACS pck file found")
    
    
    #######################################################################
    # PRIDE POS
    #######################################################################
    if type is None or type=='pride_pos':
        lGts_pride_pos = __pride_pos(ts_dir=self.dir,verbose=verbose)
        if lGts_pride_pos=={}:
            if verbose:print("-- No pride pos files found")
        else:
            lGts=lGts_pride_pos
    
    #######################################################################
    # PRIDE
    #######################################################################
    if type is None or type=='pride':
        lGts_pride = __pride(ts_dir=self.dir,verbose=verbose)
        if lGts_pride=={}:
            if verbose:print("-- No pride_files found")
        else:
            lGts=lGts_pride
    
    #######################################################################
    # MB_FILE
    #######################################################################
    if type is None or type=='mb_file':
        lGts_mb_files=__mb_files(ts_dir=self.dir,verbose=verbose)
        if lGts_mb_files=={}:
            if verbose:print("-- No mb_files found")
        else:
            lGts=lGts_mb_files
             
            try:
                self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
            except:
                print("! Warning: Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for mb_files")
     
    #######################################################################
    # TDP
    #######################################################################
    if type is None or type=='tdp':
        lGts_mb_files=__tdp(ts_dir=self.dir,verbose=verbose)
        if lGts_mb_files=={}:
            if verbose:print("-- No tdp_files found")
        else:
            lGts=lGts_mb_files
            try:
                self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
            except:
                print("! Warning: Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for mb_files")
     
    #######################################################################
    # KENV
    #######################################################################
    if type is None or type=='kenv':
        lGts_kenv=__kenv(ts_dir=self.dir,verbose=verbose)
        if lGts_kenv=={}:
            if verbose:print("-- No kenv file found")
        else:lGts=lGts_kenv
     
    #######################################################################
    # CATS
    #######################################################################
    if type is None or type=='cats':
        lGts_cats=__cats(ts_dir=self.dir,verbose=verbose)
        if lGts_cats=={}:
            if verbose:print("-- No cats file found")
        else:
            lGts=lGts_cats
            try:
                self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
            except:
                print("! Warning: Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for mb_files")
     
     
    #######################################################################
    # POS FORMAT
    #######################################################################
    if type is None or type=='pos':
        lGts_pos=__pos(ts_dir=self.dir,verbose=verbose,xyz=xyz)
        if lGts_pos=={}:
            if verbose:print("-- No Gamit/Globk pos file found")
        else:lGts=lGts_pos
    
    #######################################################################
    # TXYZ FORMAT
    #######################################################################
    if type is None or type=='txyz':
        lGts_txyz=__ts_xyz(ts_dir=self.dir,verbose=verbose)
        if lGts_txyz=={}:
            if verbose:print("-- No pyacs t_xyz file found")
        else:lGts=lGts_txyz
    
    #######################################################################
    # TRACK FORMAT
    #######################################################################
    if type is None or type=='track':
        lGts_track=__track_NEU(ts_dir=self.dir,verbose=verbose)
        if lGts_track=={}:
            if verbose:print("-- No Gamit/Globk track NEU file found")
        else:lGts=lGts_track
    
    
    for gts in list(lGts.values()):
        if gts.code not in lexclude:
            if verbose:
                print('-- adding ', gts.code)
            self.append(gts)
    
    if verbose:
        print('-- read ',len(self.lcode() ),' time series in directory ',ts_dir)
    
    
    
