"""TimePeriod.get_info implementation."""

from datetime import datetime


def get_info(self):
    str1 = datetime.isoformat(self._TimePeriod__start_time)
    str2 = datetime.isoformat(self._TimePeriod__end_time)
    s = ("%s -> %s" % (str1, str2))
    return s

