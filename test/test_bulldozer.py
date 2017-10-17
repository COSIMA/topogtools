
import os
import pytest
import subprocess as sp

class TestBulldozer():

    def test_bulldozer(self):
        my_path = os.path.abspath(os.path.dirname(__file__))
        test_data = os.path.join(my_path, 'test_data')

        bulldozer = os.path.join(my_path, '../', 'bulldozer.py')
        topog = os.path.join(test_data, 'topog.nc')
        new_topog = os.path.join(test_data, 'new_topog.nc')

        try:
            os.remove(new_topog)
        except:
            pass

        cmd = [bulldozer, topog, '--new_topog', new_topog]
        p = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        stderr, stdout = p.communicate('112, 246, 0.0, 50.0'.encode())

        if stderr != b'':
            print(stderr)
        if stdout != b'':
            print(stdout)

        assert p.returncode == 0
        assert stderr == b''
        assert stdout == b''
        assert os.path.exists(new_topog)
