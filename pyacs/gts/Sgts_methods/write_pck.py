###################################################################
def write_pck(self,outfile, verbose=True):
###################################################################
    """Write this Sgts to a pickle (.pck) file.

    Parameters
    ----------
    outfile : str
        Output file path (.pck added if no extension).
    verbose : bool, optional
        Verbose mode. Default is True.

    Returns
    -------
    None
    """
    # import
    import pickle
    import os
    from pathlib import Path



    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    outfile = Path(outfile)

    # add pck extension if not provided
     
    if outfile.suffix != '.pck':
        outfile = outfile.with_suffix('.pck')

    # create dir if it does not exist
    if os.path.dirname( outfile ) != '':
        os.makedirs( os.path.dirname( outfile ) , exist_ok=True )

    # write pickle
    ofile = open( outfile, 'wb')
    pickle.dump(self , ofile , pickle.DEFAULT_PROTOCOL)
    ofile.close()
