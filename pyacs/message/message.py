def message( ustr, level=0):
    """
    print message and log to file if configured
    """

    import logging
    from colors import red

    # Always log to file if pyacs logger is configured (regardless of level)
    pyacs_logger = logging.getLogger('pyacs')
    if pyacs_logger.handlers:
        pyacs_logger.info(ustr)

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
