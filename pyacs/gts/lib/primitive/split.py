###################################################################
def split(self, split_dates, verbose=True):
###################################################################
    """Split time series at given dates.

    Dates <= split date define the left segment. Dates are converted to
    integer seconds via at.decyear2seconds.

    Parameters
    ----------
    split_dates : list or numpy.ndarray
        Decimal-year dates for split points.
    verbose : bool, optional
        Verbose mode. Default is True.

    Returns
    -------
    list of Gts
        Gts objects from the split.
    """
    
    # import 
    import inspect
    import numpy as np
    import pyacs.gts
    import pyacs.lib.astrotime as at
    
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS ERROR
        ERROR(error)
        return []
    
    # check that data_xyz and data are consistent
    if self.data_xyz is not None:
        if self.data_xyz.shape != self.data.shape:
            ERROR(".data and .data_xyz are not consistent for time series %s. Returning empty list" % self.code)
            return []
    
    # convert split_dates to numpy array if it's a list
    if isinstance(split_dates, list):
        split_dates = np.array(split_dates)
    
    # sort split_dates to ensure proper order
    split_dates = np.sort(split_dates)
    
    # convert all dates to integer seconds
    np_date_s = at.decyear2seconds(self.data[:, 0], rounding='second')
    split_dates_s = at.decyear2seconds(split_dates, rounding='second')
    
    # create list to store split Gts objects
    split_gts_list = []
    
    # get start and end dates of the time series
    start_date_s = np_date_s[0]
    end_date_s = np_date_s[-1]
    
    # create segments based on split dates
    segment_boundaries = [start_date_s] + list(split_dates_s) + [end_date_s]
    
    if verbose:
        VERBOSE("Splitting time series %s into %d segments" % (self.code, len(segment_boundaries) - 1))
    
    # create segments
    for i in range(len(segment_boundaries) - 1):
        start_boundary = segment_boundaries[i]
        end_boundary = segment_boundaries[i + 1]
        
        # find indices for this segment
        # convention: dates <= date_offset define the left segment
        if i == 0:
            # first segment: include data <= end_boundary
            segment_indices = np.where(np_date_s <= end_boundary)[0]
        elif i == len(segment_boundaries) - 2:
            # last segment: include data > start_boundary
            segment_indices = np.where(np_date_s > start_boundary)[0]
        else:
            # middle segments: include data > start_boundary and <= end_boundary
            segment_indices = np.where((np_date_s > start_boundary) & (np_date_s <= end_boundary))[0]
        
        # check if segment has data
        if len(segment_indices) == 0:
            if verbose:
                WARNING("Segment %d (%.3f to %.3f) has no data, skipping" % 
                       (i + 1, at.datetime2decyear(at.seconds2datetime(start_boundary)), at.datetime2decyear(at.seconds2datetime(end_boundary))))
            continue
        
        # create new Gts for this segment
        new_gts = self.copy()
        
        # case .data_xyz
        if new_gts.data_xyz is not None:
            new_gts.data_xyz = self.data_xyz[segment_indices, :]
            # re-generate NEU time series
            new_gts.xyz2neu(corr=False)
        else:
            new_gts.data = self.data[segment_indices, :]
        
        # handle outliers: keep only outliers that are in this segment
        if hasattr(self, 'outliers') and self.outliers:
            segment_outliers = []
            for outlier_idx in self.outliers:
                if outlier_idx in segment_indices:
                    # adjust index for the new segment
                    new_outlier_idx = np.where(segment_indices == outlier_idx)[0][0]
                    segment_outliers.append(new_outlier_idx)
            new_gts.outliers = segment_outliers
        else:
            new_gts.outliers = []
        
        # handle offsets_dates: keep only offsets that are in this segment
        if hasattr(self, 'offsets_dates') and self.offsets_dates:
            segment_offsets = []
            for offset_date in self.offsets_dates:
                offset_date_s = at.decyear2seconds(offset_date, rounding='second')
                if start_boundary < offset_date_s <= end_boundary:
                    segment_offsets.append(offset_date)
            new_gts.offsets_dates = segment_offsets
        else:
            new_gts.offsets_dates = []
        
        # update code to indicate segment
        new_gts.code = f"{self.code}_segment_{i+1}"
        
        if verbose:
            VERBOSE("Created segment %d: %s with %d points (%.3f to %.3f)" % 
                   (i + 1, new_gts.code, len(segment_indices), 
                    at.datetime2decyear(at.seconds2datetime(start_boundary)), at.datetime2decyear(at.seconds2datetime(end_boundary))))
        
        split_gts_list.append(new_gts)
    
    if verbose:
        VERBOSE("Split completed: %d segments created" % len(split_gts_list))
    
    return split_gts_list
