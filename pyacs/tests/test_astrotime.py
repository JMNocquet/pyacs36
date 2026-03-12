"""
Unit tests for pyacs.lib.astrotime (calendar, MJD, decimal year, datetime).
"""

import numpy as np
from datetime import datetime

from pyacs.lib import astrotime as AT


# -----------------------------------------------------------------------------
# cal2mjd / mjd2cal
# -----------------------------------------------------------------------------


def test_cal2mjd_mjd2cal_round_trip():
    """Calendar and MJD round-trip (default ut=0.5)."""
    day, month, year = 29, 2, 2000
    mjd = AT.cal2mjd(day, month, year)
    d, m, y, ut = AT.mjd2cal(mjd)
    assert d == day and m == month and y == year
    assert np.isclose(ut, 0.5, atol=1e-10)


def test_cal2mjd_known_value():
    """cal2mjd(29, 2, 2000) = 51603.5 (docstring example)."""
    mjd = AT.cal2mjd(29, 2, 2000)
    assert np.isclose(mjd, 51603.5, atol=1e-6)


def test_cal2mjd_midnight():
    """cal2mjd with ut=0.0 gives integer MJD."""
    mjd = AT.cal2mjd(1, 1, 2000, ut=0.0)
    assert mjd == int(mjd)


# -----------------------------------------------------------------------------
# cal2dayno / dayno2cal
# -----------------------------------------------------------------------------


def test_cal2dayno_dayno2cal_round_trip():
    """Day of year and calendar round-trip."""
    day, month, year = 15, 6, 2020
    dayno = AT.cal2dayno(day, month, year)
    d, m = AT.dayno2cal(dayno, year)
    assert d == day and m == month


def test_cal2dayno_jan_first():
    """1 January is day of year 1."""
    assert AT.cal2dayno(1, 1, 2000) == 1
    assert AT.cal2dayno(1, 1, 2021) == 1


def test_cal2dayno_leap_day_2000():
    """29 February 2000 is day 60 (leap year)."""
    assert AT.cal2dayno(29, 2, 2000) == 60


# -----------------------------------------------------------------------------
# leap_year
# -----------------------------------------------------------------------------


def test_leap_year():
    """Leap year rules: divisible by 4, except century unless by 400."""
    assert AT.leap_year(2000)
    assert AT.leap_year(2004)
    assert not AT.leap_year(2001)
    assert not AT.leap_year(1900)
    assert AT.leap_year(2400)


# -----------------------------------------------------------------------------
# jd2mjd / mjd2jd
# -----------------------------------------------------------------------------


def test_jd2mjd_mjd2jd_round_trip():
    """JD and MJD round-trip (MJD = JD - 2400000.5)."""
    mjd = 51603.5
    jd = AT.mjd2jd(mjd)
    assert np.isclose(AT.jd2mjd(jd), mjd, atol=1e-10)


def test_jd2mjd_offset():
    """MJD is JD minus 2400000.5."""
    jd = 2451604.0
    assert np.isclose(AT.jd2mjd(jd), jd - 2400000.5, atol=1e-10)


# -----------------------------------------------------------------------------
# cal2decyear / decyear2cal
# -----------------------------------------------------------------------------


def test_cal2decyear_decyear2cal_round_trip():
    """Calendar and decimal year round-trip (midday)."""
    day, month, year = 1, 7, 2015
    decyear = AT.cal2decyear(day, month, year, ut=0.5)
    mday, m, ut = AT.decyear2cal(decyear)
    assert int(decyear) == year
    assert mday == day and m == month
    assert np.isclose(ut, 0.5, atol=1e-6)


def test_cal2decyear_start_of_year():
    """1 January midday is approximately YYYY.0."""
    decyear = AT.cal2decyear(1, 1, 2020, ut=0.5)
    assert np.isclose(decyear, 2020.0, atol=0.001)


# -----------------------------------------------------------------------------
# cal2datetime / datetime2mjd / datetime2decyear
# -----------------------------------------------------------------------------


def test_cal2datetime_midday():
    """cal2datetime at ut=0.5 gives 12:00."""
    dt = AT.cal2datetime(15, 3, 2010, ut=0.5)
    assert dt.year == 2010 and dt.month == 3 and dt.day == 15
    assert dt.hour == 12 and dt.minute == 0 and dt.second == 0


def test_cal2datetime_datetime2mjd_round_trip():
    """Calendar -> datetime -> MJD -> cal round-trip."""
    day, month, year = 20, 5, 2018
    ut = 0.25
    dt = AT.cal2datetime(day, month, year, ut=ut)
    mjd = AT.datetime2mjd(dt)
    d, m, y, u = AT.mjd2cal(mjd)
    assert d == day and m == month and y == year
    assert np.isclose(u, ut, atol=1e-6)


def test_datetime2decyear_cal2datetime_consistent():
    """datetime from cal and decyear from that datetime match cal2decyear."""
    day, month, year = 1, 6, 2012
    ut = 0.5
    decyear_from_cal = AT.cal2decyear(day, month, year, ut=ut)
    dt = AT.cal2datetime(day, month, year, ut=ut)
    decyear_from_dt = AT.datetime2decyear(dt)
    assert np.isclose(decyear_from_cal, decyear_from_dt, atol=1e-8)


# -----------------------------------------------------------------------------
# mjd2decyear / cal2decyear consistency
# -----------------------------------------------------------------------------


def test_mjd2decyear_vs_cal2decyear():
    """MJD->decyear and cal->decyear agree for same instant."""
    day, month, year = 15, 9, 2005
    ut = 0.5
    mjd = AT.cal2mjd(day, month, year, ut=ut)
    decyear_mjd = AT.mjd2decyear(mjd)
    decyear_cal = AT.cal2decyear(day, month, year, ut=ut)
    assert np.isclose(decyear_mjd, decyear_cal, atol=1e-8)
