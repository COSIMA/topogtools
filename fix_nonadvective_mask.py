#!/usr/bin/env python

import sys
import os
import numpy as np
import argparse
import netCDF4 as nc

def fix_mask(input_mask_filename, output_mask_filename):

    with nc.Dataset(input_mask_filename, 'r') as mf_in, \
         nc.Dataset(output_mask_filename, 'w') as mf_out:

        mask = mf_in.variables['mask'][:,:]
        num_lons = mf_in.dimensions['nx'].size
        num_lats = mf_in.dimensions['ny'].size
        print('mask dimensions {:11d} {:11d}'.format(num_lons, num_lats))

        # Create mask with halo (includes copied points along LH edge and RH edge and points along tripole seam)
        mask_halo = np.zeros((num_lats+1,num_lons+2), dtype=int)
        mask_halo[:-1,1:-1] = mask[:,:]

        mask_halo[:,-1] = mask_halo[:,1]
        mask_halo[:,0] = mask_halo[:,-2]
        mask_halo[-1, :] = mask_halo[-2, ::-1]

        # Find non-advective cells
        done = False
        n_total = 0
        print('non-advective cells:')
        while (not done):
            sw = (mask_halo[1:-1,  :-2] == 0) | (mask_halo[ :-2,  :-2] == 0) | (mask_halo[ :-2, 1:-1] == 0)
            se = (mask_halo[ :-2, 1:-1] == 0) | (mask_halo[ :-2, 2:  ] == 0) | (mask_halo[1:-1, 2:  ] == 0)
            ne = (mask_halo[1:-1, 2:  ] == 0) | (mask_halo[2:  , 2:  ] == 0) | (mask_halo[2:  , 1:-1] == 0)
            nw = (mask_halo[2:  , 1:-1] == 0) | (mask_halo[2:  ,  :-2] == 0) | (mask_halo[1:-1,  :-2] == 0)
            advective_cells = sw & se & ne & nw & (mask_halo[1:-1,1:-1] == 1)

            jj,ii = np.nonzero(advective_cells)
            for ij in list(zip(ii,jj+1)):
                print('{:11d} {:11d}'.format(ij[0], ij[1]))

            mask_halo[jj+1,ii+1] = 0

            n_new = np.count_nonzero(advective_cells)
            done = n_new == 0
            n_total += n_new

        print('Found {} non-advective cells'.format(n_total))

        mf_out.createDimension('nx', num_lons)
        mf_out.createDimension('ny', num_lats)

        mask_out = mf_out.createVariable('mask', 'f8', dimensions=('ny', 'nx'))
        mask_out[:,:] = mask_halo[:-1,1:-1]

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('input_mask',  help='Mask to fix.')
    parser.add_argument('output_mask', help='Output mask.')

    args = parser.parse_args()

    if not os.path.exists(args.input_mask):
        print('Error: {} was not found'.format(args.input_mask), file=sys.stderr)
        parser.print_help()
        return 1

    fix_mask(args.input_mask, args.output_mask)

if __name__ == '__main__':
    sys.exit(main())
