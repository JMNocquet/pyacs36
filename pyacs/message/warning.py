def warning( str ):
    """
    print warning message and log to file if configured
    """

    import logging
    from colors import red

    # Always log to file if pyacs logger is configured (regardless of level)
    pyacs_logger = logging.getLogger('pyacs')
    if pyacs_logger.handlers:
        pyacs_logger.warning(str)

    print(red("[PYACS WARNING] %s" % (str)))

