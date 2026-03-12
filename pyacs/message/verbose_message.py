def verbose_message( str ):
    """
    print message if global variable verbose is True and log to file if configured
    """

    import logging

    # Always log to file if pyacs logger is configured (regardless of level)
    pyacs_logger = logging.getLogger('pyacs')
    if pyacs_logger.handlers:
        pyacs_logger.info(str)
    
    # Console output controlled by level
    if logging.getLogger("my_logger").getEffectiveLevel() in [logging.INFO, logging.DEBUG]:
        print("[PYACS] %s" % (str))
