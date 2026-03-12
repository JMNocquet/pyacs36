def find_large_uncertainty( self , sigma_threshold=20 , verbose=True, lcomponent='NE', log_dir=None ):
    """Flag dates with large uncertainty as outliers.

    Parameters
    ----------
    sigma_threshold : float, optional
        Threshold in mm for flagging. Default is 20.
    verbose : bool, optional
        If True, print progress. Default is True.
    lcomponent : str or list, optional
        Components to check ('N', 'E', 'U'). Default is 'NE'.
    log_dir : str, optional
        Directory for Gts-specific log file {self.code}.log. Default is None.

    Returns
    -------
    Gts
        self (outliers updated).
    """
    
    # import
    import numpy as np
    import logging
    import os

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
        gts_logger.info(f"[{self.code}] [find_large_uncertainty] Starting method")
        gts_logger.info(f"[{self.code}] [find_large_uncertainty] Parameters: sigma_threshold={sigma_threshold}, lcomponent={lcomponent}")

    # lindex
    lidxn = []
    lidxe = []
    lidxu = []
    
    # threshold in mm
    stmm = sigma_threshold * 1.E-3
    
    # North component
    if 'N' in lcomponent:
        lidxn = np.argwhere( self.data[:,4] > stmm ).flatten().tolist()
        if gts_logger:
            gts_logger.info(f"[{self.code}] [find_large_uncertainty] North component: found {len(lidxn)} points with large uncertainty")
    
    # East component
    if 'E' in lcomponent:
        lidxe = np.argwhere( self.data[:,5] > stmm ).flatten().tolist() 
        if gts_logger:
            gts_logger.info(f"[{self.code}] [find_large_uncertainty] East component: found {len(lidxe)} points with large uncertainty")

    # Up component
    if 'U' in lcomponent:
        lidxu = np.argwhere( self.data[:,6] > stmm ).flatten().tolist()
        if gts_logger:
            gts_logger.info(f"[{self.code}] [find_large_uncertainty] Up component: found {len(lidxu)} points with large uncertainty")
    
    # update
    new_gts = self.copy()

    new_gts.outliers = list(set(sorted(self.outliers + lidxn + lidxe + lidxu)))
    
    total_outliers = len(new_gts.outliers)
    if gts_logger:
        gts_logger.info(f"[{self.code}] [find_large_uncertainty] Total outliers flagged: {total_outliers}")
        gts_logger.info(f"[{self.code}] [find_large_uncertainty] Method completed successfully")
    
    return new_gts