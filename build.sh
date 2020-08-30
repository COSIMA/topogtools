#!/bin/bash
module load intel-compiler
module load netcdf

ff="-I/$NETCDF_ROOT/include/Intel -L$NETCDF_ROOT/lib/Intel -lnetcdff -lnetcdf"
ifort check_nonadvective_mosaic.f90 -o check_nonadvective_mosaic  $ff
ifort do_partial_cells.f90  -o do_partial_cells $ff
ifort fix_nonadvective_mosaic.f90 -o fix_nonadvective_mosaic $ff
ifort float_vgrid.f90 -o float_vgrid $ff
ifort min_depth.f90 -o min_depth $ff
ifort kdtree2.f90 gen_topo.f90 -o gen_topo $ff

