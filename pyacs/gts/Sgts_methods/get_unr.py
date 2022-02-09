def get_unr(self,lcode):
    """

    :param lcode: list of GNSS site codes to be downloaded from UNR web site
    :return: new Sgts instance
    """

    # import
    from pyacs.gts.Gts import Gts
    from pyacs.gts.Sgts import Sgts

    # initialize Sgts
    ts = Sgts(read=False)


    # loop on sites

    for code in lcode:
        print("-- downloading %s from UNR" % code )
        try:
            ts.append( Gts( ).get_unr( code ) )
        except:
            print("!!!ERROR downloading %s from UNR" % code )

    # return
    return ts
