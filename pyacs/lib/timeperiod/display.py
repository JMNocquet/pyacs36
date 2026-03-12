"""TimePeriod.display implementation."""

from datetime import datetime


def display(self):
    str1 = datetime.isoformat(self._TimePeriod__start_time)
    str2 = datetime.isoformat(self._TimePeriod__end_time)
    return ("%s -> %s" % (str1, str2))

