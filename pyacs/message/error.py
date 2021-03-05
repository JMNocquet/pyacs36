def error( str , exit= False ):
    """
    print error message and optionnaly exit
    """

    from colors import red


    if exit:
        print(red("[!!!PYACS ERROR] %s. Exiting" % (str)))
        import sys
        sys.exit()

    else:
        print(red("[!!!PYACS ERROR] %s" % (str)))


