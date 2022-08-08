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

from reconstruction.modules import yardley_variants

import numpy as np
from math import radians
from scipy.spatial.transform import Rotation as R


def test_yardley_variants_nw():
    '''
    check_matrix is T1 from Table 3 of Kitahara 2004.
    '''
    ksi = [radians(5.26), radians(10.30), radians(10.5)]

    yv_test = yardley_variants(ksi)
    
    check_matrix = [
        [0.0, 0.707, -0.707],
        [-0.169, 0.697, 0.697],
        [0.986, 0.120, 0.120]
    ]

    print(f'{yv_test}, {check_matrix}')

    assert np.array_equal(yv_test[0], check_matrix)


def test_yardley_variants_ks():
    '''
    check_matrix is T1 from Table 3 of Kitahara 2006.
    '''
    ksi = [radians(0), radians(9.74), radians(9.74)]

    yv_test = yardley_variants(ksi, rot='ZYZ')
    
    check_matrix = [
        [0.742, -0.667, -0.075],
        [0.650, 0.742, -0.167],
        [0.167, 0.075, 0.983]
    ]
    print(f'{yv_test}, {check_matrix}')
    
    assert np.array_equal(yv_test[0], check_matrix)