###############################################################################
def cdata(self,data=False,data_xyz=False, tol= 0.001 ,verbose=False):
###############################################################################
    """
    Check data/data_xyz attributes

    :param data: boolean, if True, data attribute will be checked
    :param data_xyz: boolean, if True, data_xyz attribute will be checked
    :param tol: tolerance in days for two dates to be considered as the same (default 0.001 of day)
    :param verbose: boolean, verbose mode
    
    :return : boolean, True if everything is OK, False otherwise
    
    :note : in future, this routine should also whether .data and .data_xyz value are consistent
    """
    
    # import 
    import inspect
    import numpy as np
    from pyacs.gts.lib.errors import GtsInputDataNone , GtsInputData_xyzNone , GtsInputData_allNone
    from pyacs.gts.lib.errors import GtsInputDataTypeError , GtsInputDataBadDim , GtsInputDataBadNcolumns, GtsInputDataBadNrows
    from pyacs.gts.lib.errors import GtsInputData_xyzTypeError , GtsInputData_xyzBadDim , GtsInputData_xyzBadNcolumns, GtsInputData_xyzBadNrows
    from pyacs.gts.lib.errors import GtsInputDataDiffShape , GtsInputDataDiffDate 

    # CASE EITHER data OR data_xyz where asked to be checked
    
    # is there the required data
    try:
        if ( self.data is None ) and ( data ):
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( False )


    # is there the required data_xyz
    try:
        if ( self.data_xyz is None ) and ( data_xyz ):
            # raise exception
            raise GtsInputData_xyzNone(inspect.stack()[0][3],__name__,self)
    except GtsInputData_xyzNone as error:
        # print PYACS WARNING
        print( error )
        return( False )


    # GENERAL CASE

    # is there some data?
    try:
        if ( self.data is None ) and ( self.data_xyz is None ):
            # raise exception
            raise GtsInputData_allNone(inspect.stack()[0][3],__name__,self)
    except GtsInputData_allNone as error:
        # print PYACS WARNING
        print( error )
        return( False )

    # .data case
    if self.data is not None:
        
        # .data numpy array?
        try:
            if not isinstance(self.data,np.ndarray):
                # raise exception
                raise GtsInputDataTypeError(inspect.stack()[0][3],__name__,self)
        except GtsInputDataTypeError as error:
            # print PYACS WARNING
            print( error )
            return( False )

        # .data with right dimensions?
        try:
            if self.data.ndim != 2:
                # raise exception
                raise GtsInputDataBadDim(inspect.stack()[0][3],__name__,self)
        except GtsInputDataBadDim as error:
            # print PYACS WARNING
            print( error )
            return( False )
        # .data with the right shape?
        try:
            if self.data.shape[1] not in [7,10]:
                # raise exception
                raise GtsInputDataBadNcolumns(inspect.stack()[0][3],__name__,self)
        except GtsInputDataBadNcolumns as error:
            # print PYACS WARNING
            print( error )
            return( False )

        # .data with at least one observation?
        try:
            if self.data.shape[0] < 1:
                # raise exception
                raise GtsInputDataBadNrows(inspect.stack()[0][3],__name__,self)
        except GtsInputDataBadNrows as error:
            # print PYACS WARNING
            print( error )
            return( False )

    # .data_xyz case
    if self.data_xyz is not None:

        # .data numpy array?
        try:
            if not isinstance(self.data_xyz,np.ndarray):
                # raise exception
                raise GtsInputData_xyzTypeError(inspect.stack()[0][3],__name__,self)
        except GtsInputData_xyzTypeError as error:
            # print PYACS WARNING
            print( error )
            return( False )

        # .data with right dimensions?
        try:
            if self.data_xyz.ndim != 2:
                # raise exception
                raise GtsInputData_xyzBadDim(inspect.stack()[0][3],__name__,self)
        except GtsInputData_xyzBadDim as error:
            # print PYACS WARNING
            print( error )
            return( False )
        try:
            if self.data_xyz.shape[1] not in [7,10]:
                # raise exception
                raise GtsInputData_xyzBadNcolumns(inspect.stack()[0][3],__name__,self)
        except GtsInputData_xyzBadNcolumns as error:
            # print PYACS WARNING
            print( error )
            return( False )

        try:
            if self.data_xyz.shape[0] < 1:
                # raise exception
                raise GtsInputData_xyzBadNrows(inspect.stack()[0][3],__name__,self)
        except GtsInputData_xyzBadNrows as error:
            # print PYACS WARNING
            print( error )
            return( False )

    # check .data and .data_xyz consistency
    
    if ( self.data is not None ) and ( self.data_xyz is not None ):
        # check both have the same length
        try:
            if self.data.shape[0] != self.data_xyz.shape[0]:
                raise GtsInputDataDiffShape(inspect.stack()[0][3],__name__,self)
        except GtsInputDataDiffShape as error:
            # print PYACS WARNING
            print( error )
            return( False )
            
        # check both have the same dates

        try:
            if np.max((self.data[:,0]- self.data_xyz[:,0])**2) > tol / 365.25:
                raise GtsInputDataDiffDate(inspect.stack()[0][3],__name__,self)
        except GtsInputDataDiffDate as error:
            # print PYACS WARNING
            print( error )
            return( False )
            

    return(True)
