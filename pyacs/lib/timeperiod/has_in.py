"""TimePeriod.has_in implementation."""


def has_in(self, date):
    """Return whether a date falls within this period.

    Parameters
    ----------
    date : datetime
        Date to test.

    Returns
    -------
    bool
        True if date is in [start, end], False otherwise.
    """
    if not self.isdefined():
        return self.__class__()

    if self._TimePeriod__start_time <= date and self._TimePeriod__end_time >= date:
        return True
    else:
        return False

