def message( ustr, level=0):
    """
    print message
    """

    from colors import red

    if level == 0:
        print("[PYACS] %s" % ustr)

    if level == 1:
        banner = '###############################################################################'
        ustr = ustr.upper()
        print(banner)
        print("[PYACS] %s" % ustr)
        print(banner)

    if level == 2:
        banner = red('###############################################################################')
        ustr = ustr.upper()
        print(banner)
        print(banner)
        print("[PYACS] %s" % ustr)
        print(banner)
        print(banner)
