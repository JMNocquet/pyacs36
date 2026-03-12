"""
Print messages and optional logging to file.
"""

def message( message , log_file=None ):
    """
    Print a message with [PYACS] prefix and optionally write to a log file.

    Parameters
    ----------
    message : str
        Message to print.
    log_file : file-like or str, optional
        File object or path; if provided, message is appended/written.

    Returns
    -------
    None
    """


    print("[PYACS] %s " % message )

    if log_file is not None:
        import io

        if isinstance( log_file, io.IOBase):
            log_file.write("%s\n" % (message))
            return

        if isinstance( log_file, str):
            import os.path

            if os.path.isfile( log_file):
                f = open( log_file , 'a+')
                log_file.write("%s\n" % (message))
                f.close()
                return
            else:
                f = open( log_file , 'w+')
                log_file.write("%s\n" % (message))
                f.close()
                return

    return



