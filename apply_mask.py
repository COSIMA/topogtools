#!/usr/bin/env python3
"""
Apply a land mask to a topography .nc file

usage: apply_mask.py [-h] topo_in mask topo_out

positional arguments:
  topo_in     NetCDF topography input file
  mask        NetCDF land mask input file (0 on land, 1 in ocean)
  topo_out    NetCDF masked topography output file with NaN land

Andrew Kiss, July 2020
"""

try:
    import sys
    import argparse
    import xarray as xr

except ImportError:
    print('\nFatal error: modules not available.')
    print('On NCI, do the following and try again:')
    print('   module use /g/data/hh5/public/modules; module load conda/analysis3\n')
    raise

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Apply a land mask to a topography .nc file')
    parser.add_argument('topo_in', metavar='topo_in', type=str, nargs=1,
                        help='NetCDF topography input file')
    parser.add_argument('mask', metavar='mask', type=str, nargs=1,
                        help='NetCDF land mask input file (0 on land, 1 in ocean)')
    parser.add_argument('topo_out', metavar='topo_out', type=str, nargs=1,
                        help='NetCDF masked topography output file')
    args = parser.parse_args()
    topo_in = vars(args)['topo_in'][0]
    mask = vars(args)['mask'][0]
    topo_out = vars(args)['topo_out'][0]
    print('applying land mask "{}" to topography "{}"'.format(mask, topo_in))
    try:
        m = xr.open_dataset(mask).mask  # mask is 0.0 on land and 1.0 in ocean
        d = xr.open_dataset(topo_in).depth
        d = d.rename(xx='nx', yy='ny')  # to match m
        d = d.fillna(0.0)
        minimum_depth = d.where(d>0.0).min()
        d = d.where(d>0.0)  # set any zeros to nan
        d = d.fillna(minimum_depth)  # ignore any existing land mask (nan or 0)
        d = d.where(m>0.0)  # set mask to nan
        d.to_netcdf(topo_out, encoding={'depth': {'_FillValue': -9999}})
        print('saved masked topography as "{}"'.format(topo_out))
    except:
        print('Failed. Error:', sys.exc_info())
        raise

