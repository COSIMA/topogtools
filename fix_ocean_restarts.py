#!/usr/bin/env python

import sys
import os
import argparse
import shutil
import numpy as np
import netCDF4 as nc
import multiprocessing as mp
from glob import glob

"""
If the topog.nc file is changed and there are new wet or dry points then
restarts will need to be fixed up. This script helps to do this. The steps to
creating/updating restarts for the new bathymetry are as follows:

1. Do a short run from rest using the new bathymetry file.
2. Collate the above.
3. Run this script using the collated outputs from 2. as the template

Presently it only works on the 0.1 deg with 75 levels.
"""


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
            for level in range(75):
                in_data = in_fp.variables[var][0, level, :, :]
                out_data = out_fp.variables[var][0, level, :, :]
                mask = np.logical_and(np.logical_not(in_data.mask),
                                      np.logical_not(out_data.mask))
                out_data[np.where(mask)] = in_data[np.where(mask)]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("template_dir",
                        help="Name of directory containing collated restart files from a run using the new bathymetry (typically a short run from rest).")
    parser.add_argument("old_dir",
                        help="Name of directory containing collated restart files from a run using the old bathymetry.")
    parser.add_argument("output_dir",
                        help="Name of the output directory which will contain new restarts.")
    parser.add_argument("--nprocs", default=1, type=int,
                        help="Number of processes to use.")
    args = parser.parse_args()

    template_files = glob(os.path.join(args.template_dir, '*.res.nc'))
    template_files.sort()
    old_files = [os.path.join(args.old_dir, os.path.basename(f)) for f in template_files]
    output_files = [os.path.join(args.output_dir, os.path.basename(f)) for f in template_files]

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    for tf, of in zip(template_files, output_files):
        print('Copying', tf, 'to', of)
        shutil.copy(tf, of)

    pool = mp.Pool(args.nprocs)
    pool.map(copy_data, zip(old_files, output_files))

if __name__ == '__main__':
    sys.exit(main())
