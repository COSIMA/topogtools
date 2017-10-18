
import os
import numpy as np
import netCDF4 as nc
import pytest
import subprocess as sp

class TestUnmask():

    def test_unmask(self):
        my_path = os.path.abspath(os.path.dirname(__file__))
        test_data = os.path.join(my_path, 'test_data')

        unmask = os.path.join(my_path, '../', 'unmask.py')
        input = os.path.join(test_data, 'i2o.nc')
        new_input = os.path.join(test_data, 'new_i2o.nc')
        mask = os.path.join(test_data, 'kmt.nc')

        try:
            os.remove(new_input)
        except:
            pass

        cmd = [unmask, input, mask, 'kmt', '--output_file',
               new_input, '--flip_mask']
        ret = sp.call(cmd)

        assert ret == 0
