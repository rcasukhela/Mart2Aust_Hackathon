# -*- coding: utf-8 -*-
# Copyright 2018-2022 the orix developers
#
# This file is part of orix.
#
# orix is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# orix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with orix.  If not, see <http://www.gnu.org/licenses/>.

import pytest

from ....reconstruction.modules import yardley_variants

import numpy as np
from math import sin, cos
import scipy.spatial.transform.Rotation as R

@pytest.fixture
def rotation_product(ksi, string='xyz'):
    '''
    Test assumes ksi values from V.A. Yardley & E.J. Payton: 'Austenite– \
        martensite/bainite orientation relationship: characterisation \
        parameters and their application'
    '''
    
    product = R.from_euler(string,ksi[0],ksi[1],ksi[2]).as_matrix()

    return product


def test_yardley_variants():
    ksi = [5.5, 5.0, 10.5]
    
    bs_test = rotation_product(ksi)
    yv_test = modules.yardley_variants(ksi)

    assert bs_test == yv_test[1,:].reshape(3,3)