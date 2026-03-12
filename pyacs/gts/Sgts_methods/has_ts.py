###################################################################
def has_ts(self, code):
###################################################################
    """Check if a time series with the given code exists in this Sgts.

    Parameters
    ----------
    code : str
        4-character site code.

    Returns
    -------
    bool
        True if code is in Sgts, False otherwise.
    """

    if code in self.lcode():
        return True
    else:
        return False

    # def has_key(key, data):
    #     return True if key in data else False
    #
    # if has_key(code , self.__dict__):
    #     return True
    # else:
    #     return False
