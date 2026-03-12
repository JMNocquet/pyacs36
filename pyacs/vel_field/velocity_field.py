"""Velocity_Field class: composes methods from vel_field submodules."""

from pyacs.vel_field import io, info, subset, pole, common, profile, strain


class Velocity_Field:
    """Read a velocity field from a GMT psvelo file and manipulate it.

    Attributes
    ----------
    name : str, optional
        Name of the velocity field.
    file_name : str, optional
        Path to the GMT psvelo file.
    sites : list
        List of GMT_Point instances.
    """

    def __init__(self, file_name=None, name=None, lgmt_points=None, verbose=False):
        self.name = name
        self.file_name = file_name
        if lgmt_points:
            self.sites = lgmt_points
        else:
            self.sites = []


# Attach methods from submodules
Velocity_Field.read = classmethod(io.read)
Velocity_Field.write = io.write
Velocity_Field.add_point = io.add_point
Velocity_Field.remove_point = io.remove_point

Velocity_Field.info = info.info
Velocity_Field.nsites = info.nsites
Velocity_Field.l_GMT_Point = info.l_GMT_Point
Velocity_Field.print_info_site = info.print_info_site
Velocity_Field.lcode = info.lcode
Velocity_Field.site = info.site

Velocity_Field.subset = subset.subset
Velocity_Field.radial = subset.radial

Velocity_Field.calc_pole = pole.calc_pole
Velocity_Field.pole = pole.pole
Velocity_Field.substract_pole = pole.substract_pole

Velocity_Field.common = common.common

Velocity_Field.proj_profile = profile.proj_profile

Velocity_Field.strain = strain.strain
