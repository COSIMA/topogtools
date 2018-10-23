#!/usr/bin/env python

import sys
import os
import argparse
import shutil
import numpy as np
import netCDF4 as nc
import multiprocessing as mp

"""
If the topog.nc file is changed and there are new wet or dry points then
restarts will need to be fixed up. This script helps to do this. The steps to
creating/updating restarts for the new bathymetry are as follows:

1. Do a short run from rest using the new bathymetry file.
2. Collate the above.
3. Run this script

Presently it only works on the 0.1 deg.
"""

mom_restart_files = ['ocean_age.res.nc', 'ocean_bih_friction.res.nc',
                     'ocean_density.res.nc', 'ocean_sbc.res.nc',
                     'ocean_temp_salt.res.nc',
                     'ocean_velocity_advection.res.nc', 'ocean_barotropic.res.nc',
                     'ocean_con_temp.res.nc', 'ocean_frazil.res.nc',
                     'ocean_thickness.res.nc',
                     'ocean_velocity.res.nc']

def copy_data(args):

    in_f = args[0]
    out_f = args[1]
    print(in_f)
    print(out_f)

    print('Processing file {}'.format(out_f))
    with nc.Dataset(in_f) as in_fp, nc.Dataset(out_f, 'r+') as out_fp:
        for var in in_fp.variables:
            if in_fp.variables[var].shape != (1, 75, 2700, 3600):
                continue

            # Copy over level by level to save memory
            print('Processing var {}'.format(var))
            in_data = in_fp.variables[var][0, :, :, :]
            out_data = out_fp.variables[var][0, :, :, :]
            mask = np.logical_and(np.logical_not(in_data.mask),
                                  np.logical_not(out_data.mask))
            out_data[np.where(mask)] = in_data[np.where(mask)]

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("cold_start_dir",
                        help="Name of directory containing collated restart files from a cold-start run using the new bathymetry.")
    parser.add_argument("orig_dir",
                        help="Name of directory containing collated restart files from the original bathymetry.")
    parser.add_argument("output_dir",
                        help="Name of the output directory which will contain fixed restarts.")
    args = parser.parse_args()

    cold_start_files = [os.path.join(args.cold_start_dir, f) for f in mom_restart_files]
    orig_files = [os.path.join(args.orig_dir, f) for f in mom_restart_files]
    output_files = [os.path.join(args.output_dir, f) for f in mom_restart_files]

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    for csf, of in zip(cold_start_files, output_files):
        shutil.copy(csf, of)

    pool = mp.Pool(4)
    pool.map(copy_data, zip(orig_files, output_files))

if __name__ == '__main__':
    sys.exit(main())
