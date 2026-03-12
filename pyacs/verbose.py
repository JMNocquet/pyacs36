def verbose(level='INFO', print_message=True, log_dir=None):
    """Set pyacs verbose level and optionally configure file logging.

    Parameters
    ----------
    level : str, optional
        Verbose level: 'DEBUG', 'INFO', 'VERBOSE', 'WARNING', 'ERROR', or 'SILENT'.
        Default 'INFO'.
    print_message : bool, optional
        If True, print a message when the level changes. Default True.
    log_dir : str, optional
        Directory for log files. If None, file logging is not configured.

    Returns
    -------
    None
    """

    import logging
    import os
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    logging_levels={50:'CRITICAL',40:'ERROR',30:'WARNING',20:'INFO',10:'DEBUG',0:'NOTSET'}

    current_level = logging_levels[logging.getLogger("my_logger").level]

    if level.upper() not in ['DEBUG','INFO','VERBOSE','WARNING','ERROR','SILENT']:
        WARNING("level argument not understood: %s" % level)
        WARNING("level must be among: %s" % ",".join(['DEBUG','INFO','VERBOSE','WARNING','ERROR','SILENT']))
        MESSAGE("Setting verbose level to INFO")
        level='INFO'

    # Configure logging level
    if level.upper() == 'DEBUG':
        logging.getLogger("my_logger").setLevel(logging.DEBUG)

    if level.upper() == 'INFO':
        logging.getLogger("my_logger").setLevel(logging.INFO)

    if level.upper() == 'VERBOSE':
        logging.getLogger("my_logger").setLevel(logging.INFO)

    if level.upper() == 'ERROR':
        logging.getLogger("my_logger").setLevel(logging.ERROR)

    if level.upper() == 'SILENT':
        logging.getLogger("my_logger").setLevel(logging.ERROR)

    # Configure file logging if log_dir is provided
    if log_dir is not None:
        # Create log directory if it doesn't exist
        os.makedirs(log_dir, exist_ok=True)
        
        # Get the main pyacs logger
        logger = logging.getLogger('pyacs')
        
        # Remove existing file handlers to avoid duplicates
        for handler in logger.handlers[:]:
            if isinstance(handler, logging.FileHandler):
                handler.close()
                logger.removeHandler(handler)
        
        # Create log file path
        
        log_file = os.path.join(log_dir, "pyacs.log")
        
        # Configure file handler
        file_handler = logging.FileHandler(log_file, mode='a')  # append mode
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        
        # Set logger level to ensure messages are passed through
        logger.setLevel(logging.INFO)
        
        if print_message:
            MESSAGE("File logging enabled: %s" % log_file)

    if print_message:
        MESSAGE("Verbose level changed from %s to %s" % (current_level, level) )

    return