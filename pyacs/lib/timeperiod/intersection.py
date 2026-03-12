"""TimePeriod.intersection implementation."""


def intersection(self, period):
    """Return the intersection of this period with another.

    Parameters
    ----------
    period : TimePeriod
        The other time period.

    Returns
    -------
    TimePeriod
        Intersection period, or undefined if no overlap.
    """
    # If one period is not defined
    if not self.isdefined() or not period.isdefined():
        return self.__class__()

    self_start = self._TimePeriod__start_time
    self_end = self._TimePeriod__end_time

    # If period don't intersect
    if self_end < period.begin() or period.end() < self_start:
        return self.__class__()

    # For simplicity self must begin first
    if self_start > period.begin():
        return period.intersection(self)

    if self_end < period.end():
        return self.__class__(period.begin(), self_end)
    else:
        return self.__class__(period.begin(), period.end())

