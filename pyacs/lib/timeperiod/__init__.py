"""Time period representation and manipulation."""

from .errors import TimePeriodUndefinedError
from .time_period import TimePeriod

__all__ = [
    "TimePeriod",
    "TimePeriodUndefinedError",
]

