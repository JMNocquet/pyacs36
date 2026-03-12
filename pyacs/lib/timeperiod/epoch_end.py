"""TimePeriod.epoch_end implementation."""

from .errors import TimePeriodUndefinedError


def epoch_end(self):
    """Return the end time as Unix timestamp (epoch seconds).

    Returns
    -------
    int
        End time in epoch seconds.

    Raises
    ------
    TimePeriodUndefinedError
        If the period is not defined.
    """
    if not self.isdefined():
        raise TimePeriodUndefinedError("This time period is undefined")
    return int(self.end().strftime("%s"))

