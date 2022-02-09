"""
Print messages and log
"""

def message( message , log_file=None ):
    """

    :param message: message to be printed
    :param log_file: file object where message will be written
    :return:
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



