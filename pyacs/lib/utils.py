"""
Various useful routines
"""

###############################################################################
def __ensure_list_of_list(ll):
###############################################################################
    """
    Ensures ll is a list of lists
    [a,b] returns [[a,b]], and [[a,b]] returns [[a,b]]
    """

    # check this is a list
    
    if not isinstance(ll,list):
        raise TypeError('!!! __ensure_list_of_list requires a list or a list of list as argument: ',ll)
    
    # case simple list
    if not isinstance(ll[0],list):
        return([ll])
    # case list of list
    else:
        return(ll)

###############################################################################
def numpy_array_2_numpy_recarray(A,names):
###############################################################################
    """
    Converts a numpy array to a numpy recarray
    names is the names of each field
    """
    import numpy as np

    return(np.rec.array(np.core.records.array(list(tuple(A[:,:].T)),dtype={'names':names,'formats':list(map(np.dtype,A[0,:]))})))

###############################################################################
def numpy_recarray_2_numpy_array(A):
###############################################################################
    """
    Converts a structured array (with homogeneous dtype) to a np.array
    """
    import numpy as np
    return(np.array(A.view(A.dtype[0]).reshape(A.shape + (-1,))))


###################################################################
def save_np_array_with_string(A,S,fmt,outfile,comment=''):
###################################################################
    """
    saves a numpy array as a text file.
    This is equivalent to np.savetxt except that S contains strings that are added to each row.
    
    :param A: 2D numpy array to be saved
    :param S: 1D numpy string array. S.ndim = A.shape[0]
    :param fmt: format for printing
    :param outfile: out file
    
    :example:
    A = np.arange(4).reshape(-1,2)
    S = np.array(['lima','quito'])
    from pyacs.lib.utils import save_np_array_with_string
    save_np_array_with_string(A,S,"%03d %03d %s",'test.dat')
    """

    # import
    import numpy as np
    
    # check
    
    if S.size != A.shape[0]:
        raise TypeError
    # open file
    out=open(outfile,'w+')
    
    # comment
    if comment.strip() =='':
        comment = '# '
    
    if comment[0] != '#':
        comment = '# ' + comment
    
    out.write("%s\n" % comment)
    
    # loop
    
    for i in np.arange( A.shape[0] ):
        fmtstr=(fmt % tuple( list(A[i,:])+[S[i]] ))
        out.write("%s\n" %  fmtstr )
    
    # close file
    out.close()

###################################################################
def make_grid(min_lon,max_lon,min_lat,max_lat, \
              nx=None, ny=None, step_x=None, step_y=None, outfile=None, format='psvelo' , comment=''):
###################################################################
    """
    Generates a text file as a grid
    
    :param min_lon,max_lon,min_lat,max_lat: grid bounds coordinates in decimal degrees
    :param nx ny: number of points in the grid along longitude. If ny is not provided, ny=nx
    :param step_x,step_y: step for the grid. This is an alternative to nx. If step_y is not provided, step_y=step_x
    :param outfile: output file name. Default=None, no file written
    :param format: if format is None, then only lon, lat are written. If format='psvelo' (default), then the line is filled with 0. and sequentially site names
    :param comment: comment to be added to the output file.
    
    :return: the grid as 2D numpy array
    """
    
    # import
    
    import numpy as np
    from colors import red

    # check arguments
    
    if [nx,step_x] == [None,None]:
        error_str = red( "[PYACS ERROR] please provide a value for argument nx or step_x"  ) 
        print( error_str )
        raise TypeError
    
    if ( nx is not None ) and ( step_x is not None ):
        error_str = red( "[PYACS ERROR] nx or step_x argument need to be specified, not both"  ) 
        print( error_str )
        raise TypeError
    
    
    # start
    if nx is not None:
        if ny is None:
            ny = nx 
        x = np.linspace(min_lon,max_lon, nx)
        y = np.linspace(min_lat,max_lat, ny)

    if step_x is not None:
        if step_y is None:
            step_y= step_x
        x = np.arange(min_lon,max_lon, step_x)
        y = np.arange(min_lat,max_lat, step_y)

    mesh_grid = np.meshgrid(x,y)
    grid = np.vstack((mesh_grid[0].T.flatten(),mesh_grid[1].T.flatten())).T
    
    # write results
    
    if format == 'psvelo':
        output_np = np.zeros(( grid.shape[0] , 8 ) )
        output_np[:,:2] = grid
        output_np[:,7]  = np.arange( grid.shape[0] ) 

        if outfile is not None:
            np.savetxt(outfile, output_np, "%10.6lf %10.6lf %4.1lf %4.1lf %4.1lf %4.1lf %4.1lf %04d", header=comment)
        
    else:
        if outfile is not None:
            np.savetxt(outfile, grid, "%10.6lf %10.6lf", header=comment)
        
    return( grid )
    
###################################################################
def str2list_float(my_str):
###################################################################
    """
    Converts a list provided as a string to a true list of float
    
    :param my_str: string e.g. '[0,2.5,1E9]'
    :return:[0,2.5,1E9]
    """
    
    return list( map( float,  my_str.replace('[','').replace(']','').split(',') ) )   

    
