MDSTRESSLAB
===========
Installation Instructions
-------------------------

To install, do:

$ mkdir build

$ cd build

$ cmake -DCMAKE_BUILD_TYPE=Release ..

$ make

To test installation run unit tests, for example from build, do:

$ cd unit_tests/testSW

$ ./testSW

If you want to run all the unit tests, from the unit_tests directory, do:

$ ctest all  

Documentation
-------------
https://nikhil-admal.github.io/mdstresslab

Notes
-----
For Cauchy stress calculations, MDStressLab writes the stress tensor to
`<name>.stress` and the associated continuum fields to `<name>.momentum_density`
and `<name>.mass_density`.  The kinetic Cauchy stress uses atom velocities
relative to the local continuum velocity `p/rho`; Piola stress has no kinetic
contribution.  Structured grids can also be written in LAMMPS dump-grid style
with `Stress::write_voxel_grid()`, producing voxel-grid stress and, for Cauchy
stress, voxel-grid momentum and mass density files for direct OVITO
visualization.

Periodic configurations may use non-orthogonal simulation boxes.  MDStressLab
creates periodic padding atoms as needed, but it does not fold the input atoms
or grid points into the primary box.
