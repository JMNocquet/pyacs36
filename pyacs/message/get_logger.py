import logging
# function to get logger
def get_logger(log_file):
    logger = logging.getLogger('pyacs.l1trendi')
    if not logger.hasHandlers():
        logger.setLevel(logging.INFO)
        file_handler = logging.FileHandler(log_file, mode='a')  # append mode
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger
