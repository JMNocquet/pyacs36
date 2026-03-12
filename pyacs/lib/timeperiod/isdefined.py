"""TimePeriod.isdefined implementation."""


def isdefined(self):
    """Return True if both start and end times are set."""
    return self._TimePeriod__start_time is not None and self._TimePeriod__end_time is not None

