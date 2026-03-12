###################################################################
def detrend(self, method='L2', periods=[], exclude_periods=[], log_dir=None):
    ###################################################################
    """Detrend time series and store velocity (and optional offsets) in attributes.

    Parameters
    ----------
    method : str, optional
        Estimation method (e.g. 'L2'). Default is 'L2'.
    periods : list, optional
        Periods used for velocity estimate. Default is [].
    exclude_periods : list, optional
        Periods excluded from velocity estimate. Default is [].
    log_dir : str, optional
        Directory for Gts-specific log file {self.code}.log. Default is None.

    Returns
    -------
    Gts
        Detrended time series (velocity attribute set).

    Notes
    -----
    Outliers in Gts.outliers are omitted; offsets in Gts.offsets_dates are
    estimated simultaneously.
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect
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
        gts_logger.info(f"[{self.code}] [detrend] Starting method")
        gts_logger.info(f"[{self.code}] [detrend] Parameters: method={method}")
        gts_logger.info(f"[{self.code}] [detrend] Time series: {len(self.data)} points from {self.data[0, 0]:.3f} to {self.data[-1, 0]:.3f}")
        if periods:
            gts_logger.info(f"[{self.code}] [detrend] Using periods: {periods}")
        if exclude_periods:
            gts_logger.info(f"[{self.code}] [detrend] Excluding periods: {exclude_periods}")

    # after this method .data  and .data_xyz are not consistent so .data_xyz is set to None
    #self.data_xyz = None

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        if gts_logger:
            gts_logger.error(f"[{self.code}] [detrend] Input data is None")
        ERROR(error)
        return (self)
    ###########################################################################

    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend] Starting main processing")

    import copy
    outliers = copy.deepcopy(self.outliers)
    tmp_ts = self.remove_outliers()

    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend] Removed outliers, remaining points: {len(tmp_ts.data)}")

    if periods != []:
        tmp_ts = tmp_ts.extract_periods(periods)
        if gts_logger:
            gts_logger.info(f"[{self.code}] [detrend] Extracted periods, remaining points: {len(tmp_ts.data)}")
    if exclude_periods != []:
        tmp_ts = tmp_ts.exclude_periods(periods)
        if gts_logger:
            gts_logger.info(f"[{self.code}] [detrend] Excluded periods, remaining points: {len(tmp_ts.data)}")

    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend] Creating model with method={method}")

    detrended = tmp_ts.make_model(option='detrend', method=method)

    vel = detrended.velocity
    offsets_values = detrended.offsets_values

    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend] Model created successfully")
        gts_logger.info(f"[{self.code}] [detrend] Velocity estimates (N,E,U, mm/yr): [{vel[0]*1.E3:.3f}, {vel[1]*1.E3:.3f}, {vel[2]*1.E3:.3f}]")
        if len(vel) > 3:
            gts_logger.info(f"[{self.code}] [detrend] Velocity uncertainties (N,E,U, mm/yr): [{vel[3]*1.E3:.3f}, {vel[4]*1.E3:.3f}, {vel[5]*1.E3:.3f}]")
        if offsets_values is not None and len(offsets_values) > 0:
            gts_logger.info(f"[{self.code}] [detrend] Estimated {len(offsets_values)} offsets")

    new_gts = self.copy()
    new_gts.outliers = outliers
    new_gts.offsets_values = offsets_values
    new_gts.velocity = vel

    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend] Creating model time series")

    model = new_gts.mmodel()

    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend] Subtracting model from data")

    new_gts.data[:, 1:4] = new_gts.data[:, 1:4] - model.data[:, 1:4]

    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend] Method completed successfully")

    return (new_gts)
