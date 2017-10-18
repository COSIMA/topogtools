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

skip_vars = ['time']

def unmask_2d(var, mask):

    ind = nd.distance_transform_edt(mask[:, :],
                                    return_distances=False,
                                    return_indices=True)
    var[:, :] = var[tuple(ind)]


def unmask_file(filename, mask):

    def unmask_3d(v):
        for t in range(v.shape[0]):
            unmask_2d(v[t, :], mask)

    def unmask_4d(v):
        for t in range(v.shape[0]):
            unmask_3d(v[t, :])

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
                unmask_2d(var, mask)

            f.variables[v][:] = var[:]

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="""
                        The input NetCDF file with variables to unmask.""")
    parser.add_argument('mask_file', help="""
                        File containing a mask for the above variables.
                        """)
    parser.add_argument('mask_var', help="Name of mask variable.")
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

    if args.output_file is not None:
        if os.path.exists(args.output_file):
            print('Error: {} already exists'.format(args.output_file),
                  file=sys.stderr)
            parser.print_help()
            return 1

        shutil.copyfile(args.input_file, args.output_file)
        args.input_file = args.output_file

    err = unmask_file(args.input_file, mask)
    return err

if __name__ == '__main__':
    sys.exit(main())
