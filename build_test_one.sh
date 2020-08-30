#!/usr/bin/env sh
module load intel-compiler
module load netcdf
ff="-I/$NETCDF_ROOT/include/Intel -L$NETCDF_ROOT/lib/Intel -lnetcdff -lnetcdf"
ifort kdtree2.f90 test_one.f90 -o test_one $ff

