###################################################################
def write_pck(self,outfile, verbose=True):
###################################################################
    """
    writes a Sgts object as a pck (pickle)
     
    :param outfile: output file name. If not provided, a pck extension will be added.
    :param verbose: verbose mode
     
    """
    # import
    import pickle
    import os
    
    # add pck extension if not provided
     
    if outfile[-4:] != '.pck':
        outfile = outfile+'.pck'

    # create dir if it does not exist
    if os.path.dirname( outfile ) != '':
        os.makedirs( os.path.dirname( outfile ) , exist_ok=True )

    # write pickle
    ofile = open( outfile, 'wb')
    pickle.dump(self , ofile , pickle.DEFAULT_PROTOCOL)
    ofile.close()
