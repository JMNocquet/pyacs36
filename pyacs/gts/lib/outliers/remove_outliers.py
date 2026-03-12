###############################################################################
def remove_outliers(self, periods=None, in_place=False, log_dir=None):
###############################################################################
    """Remove outliers listed in self.outliers.

    Parameters
    ----------
    periods : list, optional
        If provided, restrict to these periods. Default is None.
    in_place : bool, optional
        If True, modify self; otherwise return a new Gts. Default is False.
    log_dir : str, optional
        Directory for Gts-specific log file {self.code}.log. Default is None.

    Returns
    -------
    Gts or None
        New Gts without outliers, or self if in_place=True.
    """

    # import
    import numpy as np
    import logging
    import os
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # Setup Gts-specific logging if log_dir is provided
    gts_logger = None
    if log_dir is not None:
        # Create log directory if it doesn't exist
        os.makedirs(log_dir, exist_ok=True)
        
        # Create Gts-specific log file
        log_file = os.path.join(log_dir, f"{self.code}.log")
        
        # Configure Gts-specific logger
        gts_logger = logging.getLogger(f'pyacs.gts.{self.code}')
        if not gts_logger.handlers:
            gts_logger.setLevel(logging.INFO)
            # Prevent propagation to parent loggers to avoid duplicate logging
            gts_logger.propagate = False
            file_handler = logging.FileHandler(log_file, mode='a')  # append mode
            file_handler.setLevel(logging.INFO)
            # Custom formatter with [Gts.code] prefix
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - [%(name)s] %(message)s')
            file_handler.setFormatter(formatter)
            gts_logger.addHandler(file_handler)
        
        # Log method start with [Gts.code] prefix
        gts_logger.info(f"[{self.code}] [remove_outliers] Starting method")
        gts_logger.info(f"[{self.code}] [remove_outliers] Parameters: periods={periods}, in_place={in_place}")
        gts_logger.info(f"[{self.code}] [remove_outliers] Initial outliers: {len(self.outliers)}")

    if self.outliers:
        if periods == None:
            data = np.delete(self.data, self.outliers, axis=0)
            if gts_logger:
                gts_logger.info(f"[{self.code}] [remove_outliers] Removed {len(self.outliers)} outliers from all data")
            # added JMN 2024/05/15
            if self.data_xyz is not None:
                if self.data_xyz.shape[0] == self.data.shape[0]:
                    data_xyz = np.delete(self.data_xyz, self.outliers, axis=0)
                else:
                    WARNING(".data_xyz has not the same number of rows as .data. Setting to None")
                    if gts_logger:
                        gts_logger.warning(f"[{self.code}] [remove_outliers] .data_xyz has not the same number of rows as .data. Setting to None")
                    data_xyz = None
            else:
                data_xyz = None
        else:
            lindex = self.lindex_in_periods(periods)
            ldelete = np.intersect1d(lindex, self.outliers)
            data = np.delete(self.data, ldelete, axis=0)
            if gts_logger:
                gts_logger.info(f"[{self.code}] [remove_outliers] Removed {len(ldelete)} outliers from specified periods")
            # added JMN 2024/05/15
            if self.data_xyz is not None:
                if self.data_xyz.shape[0] == self.data.shape[0]:
                    data_xyz = np.delete(self.data_xyz, ldelete, axis=0)
                else:
                    WARNING(".data_xyz has not the same number of rows as .data. Setting to None")
                    if gts_logger:
                        gts_logger.warning(f"[{self.code}] [remove_outliers] .data_xyz has not the same number of rows as .data. Setting to None")
                    data_xyz = None
            else:
                data_xyz = None

        if data.shape[0] == 0:
            ERROR("All data are outliers for site %s. Check your outlier detection threshold or remove time series." % ( self.code ))
            ERROR("Setting data as None")
            if gts_logger:
                gts_logger.error(f"[{self.code}] [remove_outliers] All data are outliers. Setting data as None")
            data = None
    else:
        data = np.copy(self.data)
        data_xyz = self.data_xyz
        if gts_logger:
            gts_logger.info(f"[{self.code}] [remove_outliers] No outliers to remove")

    new_Gts = self.copy()
    new_Gts.outliers = []
    new_Gts.data = data
    new_Gts.data_xyz = data_xyz

    if gts_logger:
        gts_logger.info(f"[{self.code}] [remove_outliers] Final data points: {data.shape[0] if data is not None else 0}")
        gts_logger.info(f"[{self.code}] [remove_outliers] Method completed successfully")

    if in_place:
        self.data = new_Gts.data.copy()
        del new_Gts
        self.outliers = []

        return (self)
    else:
        return (new_Gts)
