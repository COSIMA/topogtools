
import os
import numpy as np
import netCDF4 as nc
import pytest
import subprocess as sp

bering_straight = '112, 246, 0.0, 50.0'
red_sea = """324, 176, 0.0, 50.0
             324, 175, 0.0, 50.0
             323, 177, 0.0, 50.0
             323, 178, 0.0, 50.0
             323, 179, 0.0, 50.0
             323, 180, 0.0, 50.0
             322, 181, 0.0, 50.0
             322, 182, 0.0, 50.0
             322, 183, 0.0, 50.0
             321, 184, 0.0, 50.0
             321, 185, 0.0, 50.0
             320, 186, 0.0, 50.0
             320, 187, 0.0, 50.0
             319, 188, 0.0, 50.0
             319, 189, 0.0, 50.0
             318, 190, 0.0, 50.0
             318, 191, 0.0, 50.0
             317, 192, 0.0, 50.0
             317, 193, 0.0, 50.0
             316, 194, 0.0, 50.0"""
gulf_bay = """335, 191, 0.0, 50.0
              331, 191, 0.0, 50.0
              330, 191, 0.0, 50.0
              329, 192, 0.0, 50.0
              328, 193, 0.0, 50.0
              328, 194, 0.0, 50.0"""

class TestBulldozer():

    def test_bulldozer(self):
        my_path = os.path.abspath(os.path.dirname(__file__))
        test_data = os.path.join(my_path, 'test_data')

        bulldozer = os.path.join(my_path, '../', 'bulldozer.py')
        topog = os.path.join(test_data, 'topog.nc')
        new_topog = os.path.join(test_data, 'new_topog.nc')
        ice_mask = os.path.join(test_data, 'kmt.nc')
        new_ice_mask = os.path.join(test_data, 'new_kmt.nc')

        # Count the original number of ocean points in the ice mask.
        with nc.Dataset(ice_mask, 'r') as f:
            ocean_pts = np.sum(f.variables['kmt'][:])

        try:
            os.remove(new_topog)
        except:
            pass
        try:
            os.remove(new_ice_mask)
        except:
            pass

        cmd = [bulldozer, topog, '--new_topog', new_topog,
               '--ice_mask', new_ice_mask]
        p = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        # Fix up the Bering Straight, Red Sea and Gulf Bay
        fixes = bering_straight + '\n' + red_sea + '\n' + gulf_bay
        stderr, stdout = p.communicate(fixes.encode())
        assert stderr == b'' and stdout == b''

        assert p.returncode == 0
        assert os.path.exists(new_topog)
        assert os.path.exists(new_ice_mask)

        # Expect there to be one more ocean point than before.
        with nc.Dataset(new_ice_mask, 'r') as f:
            new_ocean_pts = np.sum(f.variables['kmt'][:])

        assert new_ocean_pts == ocean_pts + 27
