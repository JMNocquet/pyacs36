#!/usr/bin/env python
"""Launch pyacs interactive environment with common imports and objects."""

from art import tprint
import logging

import pyacs.message.message as MESSAGE


def main():
    """Print banner, set logging, import pyacs and common names. Returns a dict for interactive use."""
    tprint("PYACS")
    logging.getLogger("my_logger").setLevel(logging.WARNING)

    MESSAGE("Welcome to pyacs interactive environment")
    MESSAGE("Importing pyacs core module")
    import pyacs
    MESSAGE("pyacs version: %s" % pyacs.__version__)
    MESSAGE("Importing pyacs.gts module")
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts
    MESSAGE("Importing class Velocity_Field from pyacs.vel_field module as vf")
    from pyacs.vel_field import Velocity_Field as vf
    MESSAGE("Importing numpy as np")
    import numpy as np
    MESSAGE("Importing matplotlib.pyplot as plt")
    import matplotlib.pyplot as plt
    MESSAGE("Importing pyacs.lib.astrotime as at")
    import pyacs.lib.astrotime as at
    MESSAGE("Importing pyacs.lib.coordinates as coo")
    import pyacs.lib.coordinates as coo
    MESSAGE("Trying to read time series files")
    ts = Sgts()

    return {
        "pyacs": pyacs,
        "Gts": Gts,
        "Sgts": Sgts,
        "vf": vf,
        "np": np,
        "plt": plt,
        "at": at,
        "coo": coo,
        "ts": ts,
    }


if __name__ == "__main__":
    ns = main()
    if ns is not None:
        globals().update(ns)
