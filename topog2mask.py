#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import numpy as np
import argparse
import netCDF4 as nc

def make_mask(topog_filename, mask_filename, frac, ice=False):

    with nc.Dataset(topog_filename, 'r') as tf, \
            nc.Dataset(mask_filename, 'w') as mf:
        try:
            num_lons = tf.dimensions['xx'].size
            num_lats = tf.dimensions['yy'].size
        except KeyError:
            num_lons = tf.dimensions['nx'].size
            num_lats = tf.dimensions['ny'].size

        mf.createDimension('nx', num_lons)
        mf.createDimension('ny', num_lats)
        if ice:
            mask = mf.createVariable('kmt', 'f8', dimensions=('ny', 'nx'))
        else:
            mask = mf.createVariable('mask', 'f8', dimensions=('ny', 'nx'))
        # CICE and MOM use 0 as masked
        if (frac > 0.0):
            mask[:] = np.where((tf.variables['frac'][:] < frac) | (tf.variables['depth'][:] <= 0.0), 0, 1)
        else:
            mask[:] = np.where(tf.variables['depth'][:] <= 0.0, 0, 1)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('topog', help='The topog file.')
    parser.add_argument('cice_mask',  help='The new CICE mask file.')
    parser.add_argument('mom_mask',  help='The new MOM mask file.')
    parser.add_argument('fraction', type=float, nargs='?', default=0.0, help='Cells with a fraction of ocean lower than this value will be treated as land')
    args = parser.parse_args()

    if not os.path.exists(args.topog):
        print('Error: {} not found'.format(args.topog), file=sys.stderr)
        parser.print_help()
        return 1

    # Check that output files don't already exist
    if os.path.exists(args.cice_mask) or os.path.exists(args.mom_mask):
        print('Error: one of {} {} already exists'.format(args.cice_mask, args.mom_mask),
              file=sys.stderr)
        parser.print_help()
        return 1

    if args.fraction > 1.0 or args.fraction < 0.0:
        print('Error: fraction must be between 0 and 1, but it is {}'.format(args.fraction),
              file=sys.stderr)
        parser.print_help()
        return 1

    make_mask(args.topog, args.cice_mask, args.fraction, ice=True)
    make_mask(args.topog, args.mom_mask, args.fraction)

if __name__ == '__main__':
    sys.exit(main())
