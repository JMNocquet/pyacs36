"""Time period representation and manipulation."""

from datetime import datetime, timedelta


class TimePeriod:
    """A time period with start and end datetime.

    Attributes
    ----------
    __start_time : datetime, optional
        Start of the period.
    __end_time : datetime, optional
        End of the period.
    """

    __start_time = None
    __end_time = None

    def __init__(self, start_time=None, end_time=None):
        self.__start_time = start_time
        self.__end_time = end_time


# Attach methods (keeps original attribute privacy via name mangling).
from .isdefined import isdefined  # noqa: E402
from .begin import begin  # noqa: E402
from .epoch_begin import epoch_begin  # noqa: E402
from .end import end  # noqa: E402
from .epoch_end import epoch_end  # noqa: E402
from .intersection import intersection  # noqa: E402
from .has_in import has_in  # noqa: E402
from .display import display  # noqa: E402
from .get_info import get_info  # noqa: E402

TimePeriod.isdefined = isdefined
TimePeriod.begin = begin
TimePeriod.epoch_begin = epoch_begin
TimePeriod.end = end
TimePeriod.epoch_end = epoch_end
TimePeriod.intersection = intersection
TimePeriod.has_in = has_in
TimePeriod.display = display
TimePeriod.get_info = get_info

