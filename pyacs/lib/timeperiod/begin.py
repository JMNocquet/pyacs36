"""TimePeriod.begin implementation."""

from .errors import TimePeriodUndefinedError


def begin(self):
    """Return the start datetime of this period.

    Returns
    -------
    datetime
        Start time.

    Raises
    ------
    TimePeriodUndefinedError
        If the period is not defined.
    """
    if not self.isdefined():
        raise TimePeriodUndefinedError("This time period is undefined")
    return self._TimePeriod__start_time

