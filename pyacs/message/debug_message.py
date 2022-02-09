def debug_message( str ):
    """
    print message is log level is DEBUG
    """

    import logging

    if logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG:
        print("[PYACS] %s" % (str))
