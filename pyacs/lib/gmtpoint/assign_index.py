"""Assign index for GMT_Point."""


def assign_index(self, index):
    """Assign an index to the current GMT_Point (for linear systems).

    Parameters
    ----------
    index : int
        Index value.

    Returns
    -------
    GMT_Point
        self.
    """
    self.index = index
    return self
