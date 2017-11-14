#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import numpy as np
import argparse
import netCDF4 as nc

def make_mask(topog_filename, mask_filename, ice=False):

    with nc.Dataset(topog_filename, 'r') as tf, \
            nc.Dataset(mask_filename, 'w') as mf:
        num_lons = tf.dimensions['xx'].size
        num_lats = tf.dimensions['yy'].size
        mf.createDimension('nx', num_lons)
        mf.createDimension('ny', num_lats)
        if ice:
            mask = mf.createVariable('kmt', 'f8', dimensions=('ny', 'nx'))
        else:
            mask = mf.createVariable('mask', 'f8', dimensions=('ny', 'nx'))
        # CICE and MOM use 0 as masked
        m = np.zeros_like(mask)
        m[:] = 0
        m[np.where(tf.variables['depth'][:] > 0)] = 1
        mask[:] = m[:]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('topog', help='The topog file.')
    parser.add_argument('cice_mask',  help='The new CICE mask file.')
    parser.add_argument('mom_mask',  help='The new MOM mask file.')

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

    make_mask(args.topog, args.cice_mask, ice=True)
    make_mask(args.topog, args.mom_mask)

if __name__ == '__main__':
    sys.exit(main())
