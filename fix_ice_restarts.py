#!/usr/bin/env python

import sys
import os
import argparse
import shutil
import numpy as np
import netCDF4 as nc
import glob
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
        landmask = ~landmask

    restart = glob.glob(os.path.join(args.ice_restart_dir, 'iced.*.nc'))[0]

    skip_vars = ['time']
    unmask_file(restart, skip_vars=skip_vars)
    unmask_file(os.path.join(args.ice_restart_dir, 'i2o.nc'), missing_value=1e30, skip_vars=skip_vars)
    unmask_file(os.path.join(args.ice_restart_dir, 'o2i.nc'), missing_value=0, skip_vars=skip_vars)
    unmask_file(os.path.join(args.ice_restart_dir, 'sicemass.nc'), missing_value=1e30, skip_vars=skip_vars)
    unmask_file(os.path.join(args.ice_restart_dir, 'u_star.nc'), missing_value=1e30, skip_vars=skip_vars)

    # Apply land mask to aice. 
    apply_mask_file(restart, mask=landmask, skip_vars=skip_vars)


if __name__ == '__main__':
    sys.exit(main())
