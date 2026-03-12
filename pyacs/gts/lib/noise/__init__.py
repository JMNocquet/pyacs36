"""
Pyacs noise package: wrappers for realistic uncertainties (tsfit, CATS, add_vel_sigma).
"""

from ._noise import (
    wrms,
    realistic_sigma,
    sigma_cats,
    sigma_vel_tsfit,
    get_spectral_index,
    get_vel_and_sigma,
)
from .add_vel_sigma import add_vel_sigma

__all__ = [
    "wrms",
    "realistic_sigma",
    "sigma_cats",
    "sigma_vel_tsfit",
    "get_spectral_index",
    "get_vel_and_sigma",
    "add_vel_sigma",
]
