#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:00:47 2022

@author: paytone
"""

#def gbc_angle(q, CS, Dl, Dr, threshold, varargin):
#the inputs are: quaternion(ebsd.rotations),ebsd.CSList{p},Dl(ndx),Dr(ndx),gbcValue(p),varargin{:})



import numpy

#% get pairs of neighbouring cells {D_l,D_r} in A_D
#A_D = I_FD'*I_FD==1;
#[Dl,Dr] = find(triu(A_D,1));

# now check whether the have a misorientation heigher or lower than a threshold
m = inv(q[Dl]).*q[Dr])
criterion = np.abs(dot(m,quaternion.id)) > np.cos(threshold/2.0);
if any(~criterion):
  qcs = quaternion(CS.properGroup)
  criterion(~criterion) = np.amax(np.abs(dot_outer(m(~criterion),qcs)),[],2) > np.cos(threshold/2.);

# o_Dl = orientation(q(Dl),CS,symmetry);
# o_Dr = orientation(q(Dr),CS,symmetry);
# criterion = dot(o_Dl,o_Dr) > cos(threshold/2);

return criterion
