###################################################################
def correct_duplicated_dates(self,action='correct',tol= .1, in_place=False,verbose=False):
###################################################################
    """
    Check or remove duplicated dates in a time series
    
    :param action: 'correct' (default) or 'check'
    :param tol: tolerance for two dates to be considered as the same (default = 0.1 day)
    :param in_place: boolean, if True, 
    :param verbose: verbose mode
    """

    # import
    import inspect
    import numpy as np
    import pyacs.gts.Gts
    from pyacs.gts.lib.errors import GtsCDataError

    # check data
    try:
        if not self.cdata(verbose=verbose):
            raise GtsCDataError(inspect.stack()[0][3],__name__,self)
    except GtsCDataError as error:
        print( error )
        return( self )

    # new_gts
    
    new_gts=self.copy( )

    # reorder first
    
    self.reorder(verbose=verbose)
    
    # check duplicate
    
    if self.data is not None:
        diff=np.diff(self.data[:,0])
        lindex = np.where(diff < tol / 365.25 )

        # there are duplicated dates
        if lindex[0].shape[0] > 0:
            if verbose:
                print('-- There are ',len(lindex),' duplicated dates for .data at:')
                
            # some dates are duplicated

            new_data=np.copy(self.data)
            new_data_xyz=np.copy(self.data_xyz)
            
            if action == 'check':
                if not verbose:
                    print('-- There are duplicated dates for .data at:')
                for index in lindex:
                    print('-- ',self.data[index,0],' --- ',self.data[index+1,0])
                
                
            if action == 'correct':
                if verbose:
                    print('-- deleting ',lindex[0].shape[0],' entries in Gts')
                # keep the second date
                new_data=np.delete(self.data, lindex, axis=0)
                if new_data_xyz is not None:
                    new_data_xyz=np.delete(self.data_xyz, lindex, axis=0)
        
        # there are no duplicated dates
        else:
            new_data=np.copy(self.data)
            new_data_xyz=np.copy(self.data_xyz)
            if verbose:
                print('-- no duplicated dates in .data for ',self.code)
            
    
    # deal with outliers. assume that they are defined for .data

    if ( lindex[0].shape[0] > 0 ) and ( action == 'correct' ):
        
        if verbose:
            print('-- updating outliers list for ',self.code)
        
        ldate_outliers=self.data[:,0][self.outliers]
        ldate_removed=self.data[:,0][lindex]
    
        ldate_tmp=list( set(ldate_outliers).intersection(set(ldate_removed)) )
        
        if verbose:
            print('-- removing ',len(ldate_tmp),' outliers from Gts.outliers')
        
        lupdated_outliers=pyacs.gts.Gts.get_index_from_dates(ldate_tmp, new_data, tol=0.1)
    
    else:
        lupdated_outliers=self.outliers
        
    # in_place ?
    
    if in_place:
        if self.data is not None:
            self.data=new_data
            # update outliers
            self.outliers=lupdated_outliers
            
        if self.data_xyz is not None:self.data_xyz=new_data_xyz
        return(self)

    else:
        if new_gts.data is not None:new_gts.data=new_data
        if new_gts.data_xyz is not None:new_gts.data_xyz=new_data_xyz
        new_gts.outliers=lupdated_outliers
        return(new_gts)
