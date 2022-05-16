#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:00:47 2022

@author: paytone
"""

def householderMatrix(v):
    import numpy as np
    v = np.atleast_2d(v).T
    H = np.eye(3) - 2. / np.matmul(v.T, v) * np.matmul(v, v.T)
    return H

def translationMatrix(s):
    import numpy as np
    T  = np.array([[ 1., 0., s[0]],[0., 1., s[1]],[0., 0., 1.]])
    return T


import numpy as np


T_out = np.loadtxt('translation_output.txt', delimiter=',')
T_in = np.loadtxt('translation_input.txt', delimiter=',')

T_test = translationMatrix(T_in)

H_out = np.loadtxt('householder_output.txt', delimiter=',')
H_in = np.loadtxt('householder_input.txt', delimiter=',')

H_test = householderMatrix(H_in)