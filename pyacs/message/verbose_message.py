def verbose_message( str ):
    """
    print message is global variable verbose is True
    """

    import logging

    if logging.getLogger("my_logger").getEffectiveLevel() in [logging.INFO, logging.DEBUG]:
        print("[PYACS] %s" % (str))
