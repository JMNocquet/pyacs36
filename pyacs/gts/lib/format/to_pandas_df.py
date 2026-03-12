"""
Convert pyacs Gts instance to pandas DataFrame.
"""

def to_pandas_df( self , data_xyz=False , uncertainty=False , round=False ):
    """
    Convert a pyacs Gts to a pandas DataFrame.

    Parameters
    ----------
    data_xyz : bool, optional
        If True, use .data_xyz (X,Y,Z); otherwise use .data (N,E,U).
    uncertainty : bool, optional
        If True, include uncertainty columns.
    round : bool, optional
        If True, round index (unused).

    Returns
    -------
    pandas.DataFrame
        Index is datetime; columns are position and optionally uncertainties.

    Notes
    -----
    Uncertainties are included only when uncertainty=True.
    """

    import pandas as pd
    import pyacs.lib.astrotime as at

    if data_xyz:
        if uncertainty:
            df = pd.DataFrame(self.data_xyz[:, 1:10], columns=['X', 'Y', 'Z','SX','SY','SZ','SXY','SXZ','SYZ'])
        else:
            df = pd.DataFrame(self.data_xyz[:, 1:4], columns=['X', 'Y', 'Z'])

    else:
        if uncertainty:
            df = pd.DataFrame(self.data[:, 1:10], columns=['N', 'E', 'U','SN','SE','SU','SNE','SNU','SEU'])
        else:
            df = pd.DataFrame(self.data[:, 1:4], columns=['N', 'E', 'U'])

    df.index = pd.to_datetime(at.decyear2datetime(self.data[:, 0]))

    return df
