#!/usr/bin/env python

import sys
import shutil
import git
import argparse
import datetime
import numpy as np
import netCDF4 as nc

def create_history_record(repo_dir=None):
    """
    Create a new history record that documents how this command was run.
    """

    if not repo_dir:
        repo_dir = os.path.dirname(os.path.realpath(__file__))

    time_stamp = datetime.datetime.now().isoformat()
    exe = sys.executable
    args = " ".join(sys.argv)
    repo = git.Repo(repo_dir)
    git_url = str(list(repo.remote().urls)[0])
    git_hash = str(repo.heads[0].commit)[0:7]

    entry = "{}: {} {} (Git URL: {}, Git hash: {})".format(time_stamp, exe, args, git_url, git_hash)

    return entry


def get_neighbours(pt, max_j, max_i, exclude):

    j, i = pt

    ngbs = [(j, i+1), (j+1, i+1),
            (j+1, i), (j+1, i-1),
            (j, i-1), (j-1, i-1),
            (j-1, i), (j-1, i+1)]

    def within_bounds(p):
        if p[0] >= max_j or p[0] < 0 or p[1] >= max_i or p[1] < 0:
            return False
        else:
            return True

    ngbs = list(filter(within_bounds, ngbs))
    ngbs = list(filter(lambda n: n not in exclude, ngbs))

    return ngbs

def find_connected_points(depth, sj, si, cut_or_fill_depth, fill=False):

    connected = set([(sj, si)])
    visited = set([(sj, si)])

    max_j = depth.shape[0]
    max_i = depth.shape[1]
    to_visit = get_neighbours((sj, si), max_j, max_i, [])

    while len(to_visit) != 0:
        new_to_visit = set([])

        for p in to_visit:
            visited.add(p)
            if fill:
                if depth[p] >= cut_or_fill_depth:
                    connected.add(p)
                    new_to_visit.update(get_neighbours(p, max_j, max_i, visited))
            else:
                if depth[p] <= cut_or_fill_depth:
                    connected.add(p)
                    new_to_visit.update(get_neighbours(p, max_j, max_i, visited))

        to_visit = new_to_visit

    return connected

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("start_i", help="i coordinate of a point in basin/island.", type=int)
    parser.add_argument("start_j", help="j coordinate of a point in basin/island.", type=int)
    parser.add_argument("input_file", help="The input bathymetry file.")
    parser.add_argument("output_file", help="The output bathymetry file.")
    parser.add_argument("--cut_or_fill_depth", default=None, type=float,
                        help="The depth to cut (or fill) to.")
    parser.add_argument("--fill", action='store_true', default=False,
                        help="Fill instead of cutting.")
    args = parser.parse_args()

    shutil.copy(args.input_file, args.output_file)

    f = nc.Dataset(args.output_file, mode='r+')
    depth = f.variables['depth'][:]

    if not args.cut_or_fill_depth:
        args.cut_or_fill_depth = depth[args.start_j, args.start_i]

    connected_points = find_connected_points(depth, args.start_j, args.start_i,
                                             args.cut_or_fill_depth, args.fill)

    depth[list(zip(*connected_points))] = args.cut_or_fill_depth

    f.variables['depth'][:] = depth

    history = create_history_record()
    if hasattr(f, 'history'):
        f.history = '{} \n {}'.format(f.history, history)
    else:
        f.history = history

    f.close()


if __name__ == '__main__':
    sys.exit(main())
