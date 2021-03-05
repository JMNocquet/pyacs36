def message( str , level=0 ):
    """
    print message
    """

    if level > 0:
        banner = '###############################################################################'
        str = str.upper()
        print(banner)
        print("[PYACS] %s" % (str))
        print(banner)

    else:
        print("[PYACS] %s" % (str))

