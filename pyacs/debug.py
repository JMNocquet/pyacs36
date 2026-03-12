def debug():
    """

    :return: True if verbose mode is debug, returns False otherwise
    """

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    logging_levels={50:'CRITICAL',40:'ERROR',30:'WARNING',20:'INFO',10:'DEBUG',0:'NOTSET'}

    if logging_levels[logging.getLogger("my_logger").level]=='DEBUG':
        return True
    else:
        return False
