MDSTRESSLAB
===========
Installation Instructions

To install, from mdstresslab/ do:

$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make

To test installation run unit tests, for example from build/ do:

$ cd unit_tests/testSW
$ ./testSW

If you want to run all the unit tests, from the unit_tests directory, do:

$ ctest all  
