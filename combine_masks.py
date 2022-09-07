#!/usr/bin/env python

import sys
import os
import numpy as np
import argparse
import netCDF4 as nc

def fix_mask(mask1_filename, mask2_filename, grid_filename, lat, output_mask_filename):

    with nc.Dataset(mask1_filename, 'r') as mf1, \
         nc.Dataset(mask2_filename, 'r') as mf2, \
         nc.Dataset(grid_filename, 'r') as gr, \
         nc.Dataset(output_mask_filename, 'w') as mf_out:

        # Sanity checks
        num_lons = mf1.dimensions['nx'].size
        num_lats = mf1.dimensions['ny'].size
        if 2*num_lons != gr.dimensions['nx'].size or 2*num_lats != gr.dimensions['ny'].size:
            print('Error: dimensions of {} and {} are not compatible'.format(mask1_filename, grid_filename), file=sys.stderr)
            return 1

        if num_lons != mf2.dimensions['nx'].size or num_lats != mf2.dimensions['ny'].size:
            print('Error: dimensions of {} and {} are not compatible'.format(mask1_filename, mask2_filename), file=sys.stderr)
            return 1

        # Get the coordinates of the T points from the supergrid:
        yt = gr.variables['y'][1::2,1::2]

        # Combine the two masks, taking mask 1 north of 'lat' and mask 2 south of 'lat'.
        combined = np.where(yt > lat, mf1.variables['mask'][:], mf2.variables['mask'][:])

        # Output mask
        mf_out.createDimension('nx', num_lons)
        mf_out.createDimension('ny', num_lats)
        mask_out = mf_out.createVariable('mask', 'f8', dimensions=('ny', 'nx'))
        mask_out[:] = combined[:]

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('mask1', help='First mask to combine.')
    parser.add_argument('mask2', help='Second mask to combine.')
    parser.add_argument('grid', help='Ocean super-grid.')
    parser.add_argument('output_mask', help='The combined mask.')
    parser.add_argument('latitude', type=float, help='Take mask1 north of this value and mask2 otherwise.')

    args = parser.parse_args()

    if not os.path.exists(args.mask1) or not os.path.exists(args.mask2) or not os.path.exists(args.grid):
        print('Error: one of {}, {} or {} was not found'.format(args.topog, args.input_mask, args.grid), file=sys.stderr)
        parser.print_help()
        return 1

    return fix_mask(args.mask1, args.mask2, args.grid, args.latitude, args.output_mask)

if __name__ == '__main__':
    sys.exit(main())
