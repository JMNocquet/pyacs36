"""
Converts a pyacs gts instance to a pandas DataFrame
"""

def to_pandas_df(self):
    """
    Converts a pyacs Gts to a pandas dataframe

    :return: pandas DataFrame

    :note: uncertainties are not imported.
    """

    import pandas as pd
    import pyacs.lib.astrotime as at


    df = pd.DataFrame(self.data[:, 1:4], columns=['N', 'E', 'U'])
    df.index = pd.to_datetime(at.decyear2datetime(self.data[:, 0]))

    return df
