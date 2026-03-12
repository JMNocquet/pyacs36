"""TimePeriod.end implementation."""

from .errors import TimePeriodUndefinedError


def end(self):
    """Return the end datetime of this period.

    Returns
    -------
    datetime
        End time.

    Raises
    ------
    TimePeriodUndefinedError
        If the period is not defined.
    """
    if not self.isdefined():
        raise TimePeriodUndefinedError("This time period is undefined")
    return self._TimePeriod__end_time

