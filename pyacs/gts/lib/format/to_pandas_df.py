"""
Converts a pyacs gts instance to a pandas DataFrame
"""

def to_pandas_df( self , data_xyz=False , uncertainty=False , round=False ):
    """
    Converts a pyacs Gts to a pandas dataframe

    :return: pandas DataFrame

    :note: uncertainties are not imported.
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
