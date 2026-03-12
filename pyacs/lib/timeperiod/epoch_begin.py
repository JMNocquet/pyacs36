"""TimePeriod.epoch_begin implementation."""

from .errors import TimePeriodUndefinedError


def epoch_begin(self):
    """Return the start time as Unix timestamp (epoch seconds).

    Returns
    -------
    int
        Start time in epoch seconds.

    Raises
    ------
    TimePeriodUndefinedError
        If the period is not defined.
    """
    if not self.isdefined():
        raise TimePeriodUndefinedError("This time period is undefined")
    return int(self.begin().strftime("%s"))

