def error( str , exit= False ):
    """
    print error message and optionally exit, also log to file if configured
    """

    import logging
    from colors import red

    # Always log to file if pyacs logger is configured (regardless of level)
    pyacs_logger = logging.getLogger('pyacs')
    if pyacs_logger.handlers:
        if exit:
            pyacs_logger.error(str + " - Exiting")
        else:
            pyacs_logger.error(str)

    if exit:
        print(red("[!!!PYACS ERROR] %s. Exiting" % (str)))
        import sys
        sys.exit()

    else:
        print(red("[!!!PYACS ERROR] %s" % (str)))


