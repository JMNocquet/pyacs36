Velocity field
==============

The **pyacs.vel_field** package provides the :class:`~pyacs.vel_field.Velocity_Field` class to read, manipulate, and write velocity fields in GMT psvelo format. Each site is represented as a :class:`~pyacs.lib.gmtpoint.GMT_Point` with position, velocity components (e.g. Ve, Vn), and uncertainties.

**Full API reference:** :doc:`api/pyacs.vel_field`

**Backward compatibility:** :class:`~pyacs.vel_field.Velocity_Field` is also available as :class:`~pyacs.lib.vel_field.Velocity_Field` (``pyacs.lib.vel_field`` is a thin wrapper that imports from ``pyacs.vel_field``).


Velocity_Field class API
------------------------

Full reference: :class:`pyacs.vel_field.Velocity_Field`

The **Velocity_Field** class holds a collection of GMT_Point sites (lon, lat, Ve, Vn, uncertainties) and provides methods to read/write GMT psvelo files, subset sites, compute Euler poles, remove pole predictions, project along profiles, and estimate strain rates.


Summary of available methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

I/O
^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.vel_field.Velocity_Field.read`
     - Read a velocity field from a GMT psvelo file (class method).
   * - :meth:`~pyacs.vel_field.Velocity_Field.write`
     - Write the velocity field to a GMT psvelo file.
   * - :meth:`~pyacs.vel_field.Velocity_Field.add_point`
     - Append a GMT_Point to the field.
   * - :meth:`~pyacs.vel_field.Velocity_Field.remove_point`
     - Remove a site by code.


Info and accessors
^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.vel_field.Velocity_Field.info`
     - Print basic information about the velocity field.
   * - :meth:`~pyacs.vel_field.Velocity_Field.print_info_site`
     - Print information for one site by code.
   * - :meth:`~pyacs.vel_field.Velocity_Field.nsites`
     - Return the number of sites.
   * - :meth:`~pyacs.vel_field.Velocity_Field.lcode`
     - Return list of site codes.
   * - :meth:`~pyacs.vel_field.Velocity_Field.site`
     - Return a GMT_Point by code.
   * - :meth:`~pyacs.vel_field.Velocity_Field.l_GMT_Point`
     - Return the field as a list of GMT_Point instances.


Subset and transform
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.vel_field.Velocity_Field.subset`
     - Return a new velocity field with a subset of sites (lonly / lexclude).
   * - :meth:`~pyacs.vel_field.Velocity_Field.radial`
     - Return a field with radial and tangential components about a center.


Euler pole
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.vel_field.Velocity_Field.pole`
     - Compute Euler pole from the velocity field (WLS or L1).
   * - :meth:`~pyacs.vel_field.Velocity_Field.calc_pole`
     - Euler pole calculation for multiple plates; writes euler_sum.dat and per-plate outputs.
   * - :meth:`~pyacs.vel_field.Velocity_Field.substract_pole`
     - Subtract (or add) the prediction of an Euler pole from the field.


Common sites and profile
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.vel_field.Velocity_Field.common`
     - Return sites common to this field and a list of GMT_Points (e.g. from SINEX).
   * - :meth:`~pyacs.vel_field.Velocity_Field.proj_profile`
     - Project velocity components along a great-circle profile.


Strain
^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Method
     - Description
   * - :meth:`~pyacs.vel_field.Velocity_Field.strain`
     - Calculate strain rate from a list of site codes (velocity gradient tensor, principal axes, rotation).


API reference
~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   :caption: Velocity field API

   api/pyacs.vel_field
