#!/usr/bin/env python

import sys
import os
import argparse
import shutil
import numpy as np
import netCDF4 as nc
import multiprocessing as mp

from unmask import unmask_file, apply_mask_file

"""
If the ice layout / block size is changed then things can stop working (FIXME: I don't really understand why).

CICE uses the ice fraction field to decided where to run.

This script fixes up the restart fields so that there is no ice over land.
"""

cice_restart_files = ['i2o.nc', 'iced.*.nc', 'monthly_sstsss.nc',  'o2i.nc',  'sicemass.nc',  'u_star.nc']


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("ice_restart_dir",
                        help="Name of directory containing ice restart files.")
    args = parser.parse_args()

    with nc.Dataset(os.path.join(args.ice_restart_dir, 'kmt.nc')) as f:
        landmask = np.array(f.variables['kmt'][:], dtype=bool)

    restart = glob.glob(os.path.join(args.ice_restart_dir, 'iced.*.nc'))[0]

    # FIXME: use a missing value instead of a mask for this.
    umask_file(restart, landmask, skip_vars=) 
    umask_file(os.path.join(args.ice_restart_dir, 'i2o.nc'), landmask) 
    umask_file(os.path.join(args.ice_restart_dir, 'monthly_sstsss.nc'), landmask) 
    umask_file(os.path.join(args.ice_restart_dir, 'o2i.nc'), landmask) 
    umask_file(os.path.join(args.ice_restart_dir, 'sicemass.nc'), landmask) 
    umask_file(os.path.join(args.ice_restart_dir, 'u_star.nc'), landmask) 

    # Apply mask to aice. 
    apply_mask_file(restart, 


if __name__ == '__main__':
    sys.exit(main())
