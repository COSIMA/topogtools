#!/usr/bin/env python

import sys
import os
import time
import shutil
import numpy as np
import argparse
import netCDF4 as nc
from scipy import ndimage as nd

"""
Remove masks from all fields in a file using nearest neighbour.
"""


def unmask_2d(var, mask, missing_value):

    if mask is None:
        assert missing_value is not None, 'Need mask or missing_value'
        mask = np.zeros_like(var)
        mask[np.where(var >= missing_value)] = 1

    ind = nd.distance_transform_edt(mask[:, :],
                                    return_distances=False,
                                    return_indices=True)
    var[:, :] = var[tuple(ind)]
    print('2d done', flush=True)


def unmask_3d(v, mask, missing_value):
    for t in range(v.shape[0]):
        unmask_2d(v[t, :], mask, missing_value)

def unmask_4d(v, mask, missing_value):
    for t in range(v.shape[0], mask, missing_value):
        unmask_3d(v[t, :])


def unmask_file(filename, mask=None, missing_value=None, skip_vars=[]):

    with nc.Dataset(filename, 'r+') as f:
        for v in f.variables:
            if v in skip_vars:
                continue

            var = f.variables[v][:]

            if len(var.shape) == 4:
                unmask_4d(var)
            elif len(var.shape) == 3:
                unmask_3d(var)
            else:
                unmask_2d(var, mask, missing_value)

            f.variables[v][:] = var[:]

    return 0


def apply_mask_2d(v, landmask, mask_val):
    v[np.where(landmask)] = mask_val


def apply_mask_3d(v, landmask, mask_val):
    for d in range(v.shape[0]):
        apply_landmask_2d([v[d, :], landmask, mask_val)

def apply_mask_4d(v, landmask, mask_val):
    for t in range(v.shape[0], landmask, mask_val):
        apply_landmask_3d(v[t, :])


def apply_mask_file(filename, mask, mask_val=0.0, skip_vars=[]):

    with nc.Dataset(filename, 'r+') as f:
        for v in f.variables:
            if v in skip_vars:
                continue

            var = f.variables[v][:]

            if len(var.shape) == 4:
                unmask_4d(var)
            elif len(var.shape) == 3:
                unmask_3d(var)
            else:
                unmask_2d(var, mask, missing_value)

            f.variables[v][:] = var[:]

    return 0




def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="""
                        The input NetCDF file with variables to unmask.""")
    parser.add_argument('--mask_file', default=None, type=str, help="""
                        File containing a mask for the above variables.
                        """)
    parser.add_argument('--mask_var', default=None, type=str, help="Name of mask variable.")
    parser.add_argument('--missing_value', default=None, type=float, help="Missing value to use as mask.")
    parser.add_argument('--output_file', default=None, help="""
                        Name of the output file. If not given the input will
                        be modified in place""")
    parser.add_argument('--flip_mask', action='store_true', default=False,
                        help="""Flip the mask convenction.
                                By default 1 is land""")

    args = parser.parse_args()

    # Check that files exist
    if not os.path.exists(args.input_file):
        print('Error: {} not found'.format(args.input_file),
              file=sys.stderr)
        parser.print_help()
        return 1

    if args.mask_file is not None:
        assert args.mask_var is not None, "Please supply --mask_var"
        if not os.path.exists(args.mask_file):
            print('Error: {} not found'.format(args.mask_file),
                  file=sys.stderr)
            parser.print_help()
            return 1
        else:
            with nc.Dataset(args.mask_file) as f:
                if args.mask_var not in f.variables:
                    print('Error: var {} not found in {}'.format(args.mask_var,
                           args.mask_file), file=sys.stderr)
                    parser.print_help()
                    return 1
                else:
                    mask = f.variables[args.mask_var][:]
                    if args.flip_mask:
                        mask[:] = abs(mask[:] - 1)
    else:
        assert args.missing_value is not None, "Please supply --mask_file or --missing_value"
        assert args.mask_file is None
        assert args.mask_var is None
        mask = None


    if args.output_file is not None:
        if os.path.exists(args.output_file):
            print('Error: {} already exists'.format(args.output_file),
                  file=sys.stderr)
            parser.print_help()
            return 1

        shutil.copyfile(args.input_file, args.output_file)
        args.input_file = args.output_file

    skip_vars = ['time', 'Time', 'zaxis_1', 'xaxis_1', 'yaxis_1', 'zaxis_2']
    err = unmask_file(args.input_file, mask, args.missing_value, skip_vars=skip_vars)
    return err

if __name__ == '__main__':
    sys.exit(main())
