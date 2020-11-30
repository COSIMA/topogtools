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
If a change to the processor layout alters the cpu mask then
restarts will need to be fixed up. This script helps to do this. The steps to
creating/updating restarts for the new layout are as follows:

1. Do a short run from rest using the new layout file.
2. Collate the above.
3. Run this script using the collated outputs from 2. as the template
"""


def copy_mask(args):

    in_f, out_f = args

    print('\nProcessing file {}'.format(out_f))
    print(  '     using file {}'.format(in_f))
    with nc.Dataset(in_f) as in_fp, nc.Dataset(out_f, 'r+') as out_fp:
        for var in in_fp.variables:
            varshape = in_fp.variables[var].shape
            if len(varshape) < 2:
                print('     --- skipping var {} with shape {}'.format(var, varshape))
            elif len(varshape) < 4:
                print('     Processing var {} with shape {}'.format(var, varshape))
                in_data = in_fp.variables[var][...]
                out_data = out_fp.variables[var][...]
                unmasked = np.logical_and(out_data.mask, np.logical_not(in_data.mask))  # newly unmasked cells
                out_data.mask = in_data.mask
                out_data[np.where(unmasked)] = in_data[np.where(unmasked)]
                out_fp.variables[var][...] = out_data
            else:
                # Copy over level by level to save memory.
                nlevels = varshape[-3]  # Assume 3rd-last index is level.
                print('     Processing var {} with shape {} by looping over {} levels'.format(var, varshape, nlevels))
                for level in range(nlevels):
                    in_data = in_fp.variables[var][..., level, :, :]
                    out_data = out_fp.variables[var][..., level, :, :]
                    unmasked = np.logical_and(out_data.mask, np.logical_not(in_data.mask))  # newly unmasked cells
                    out_data.mask = in_data.mask
                    out_data[np.where(unmasked)] = in_data[np.where(unmasked)]
                    out_fp.variables[var][..., level, :, :] = out_data


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("template_dir",
                        help="Name of directory containing collated restart files from a run using the new layout (typically a short run from rest).")
    parser.add_argument("old_dir",
                        help="Name of directory containing collated restart files from a run using the old layout.")
    parser.add_argument("output_dir",
                        help="Name of the output directory which will contain new restarts.")
    parser.add_argument("--nprocs", default=1, type=int,
                        help="Number of processes to use.")
    args = parser.parse_args()

    template_files = glob(os.path.join(args.template_dir, '*.res.nc'))
    template_files.sort()
    old_files = [os.path.join(args.old_dir, os.path.basename(f)) for f in template_files]
    output_files = [os.path.join(args.output_dir, os.path.basename(f)) for f in template_files]

    print('\ntemplate_files =', template_files)
    print('\nold_files =', old_files)
    print('\noutput_files =', output_files)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print()
    for oldf, outf in zip(old_files, output_files):
        print('Copying', oldf, 'to', outf)
        shutil.copy2(oldf, outf)

    pool = mp.Pool(args.nprocs)
    pool.map(copy_mask, zip(template_files, output_files))


if __name__ == '__main__':
    sys.exit(main())
