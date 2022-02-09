###################################################################
def sel_period(self,period,min_data=2,verbose=True):
###################################################################
    """
    selects time series having some data for a given period
     
    :param period: [start,end], start and end period as decimal years
    :param min_data: minimum number of data for a time series to be kept
    :pram verbose: verbose mode
     
    :return: a new Sgts instance
    """
    # debug
    debug = False

    # loop on gts
    lsite = []
    for code in sorted( self.lcode() ):
        if debug:
            print("-- Testing %s" % code)
        wts = self.__dict__[code].extract_periods( period )
        if wts.data is not None:
            if wts.data.shape[0] >= min_data:
                lsite.append(code)
                if debug: print("-- Keeping GPS sites:%s " % code )
            else:
                if verbose:
                    print("-- Less than %d data in requested period: %s" % ( min_data , code ) )
        else:
            if verbose:
                print("-- No data for requested period site:%s" %  code )

    if verbose:
        print("-- %d over %d sites kept" % ( len(lsite) , len(self.lcode()) ))            
    return self.sub(linclude = lsite )            
