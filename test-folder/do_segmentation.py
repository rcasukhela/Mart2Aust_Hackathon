#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Tyler Martin
Github: martint98
Date: 5/16/2022
"""

import numpy as np
from scipy import sparse
from spatialDecompFunctions import gbc_angle


def do_segmentation(I_FD, ebsd, varargin):
    # Output
    # A_Db - adjecency matrix of grain boundaries
    # A_Do - adjecency matrix inside grain connections

    ## if numel(gbcValue) == 1 && length(ebsd.CSList) > 1
    if np.size(gbcValue) == 1 and max((ebsd.CSList).shape) > 1:
        ##   gbcValue = repmat(gbcValue,size(ebsd.CSList))
        gbcValue = np.repeat(gbcValue, np.shape(ebsd.CSList))

    # get pairs of neighbouring cells {D_l,D_r} in A_D
    ## A_D = I_FD'*I_FD==1
    A_D = np.atleast_2d(I_FD).T.conj()
    ## [Dl,Dr] = find(triu(A_D,1))
    Dl, Dr = np.triu(A_D, 1)          # Get upper triangular part of matrix
    # exampleArray > 50  # Example 0fbBoolean masking using greater-than comparator if even needed

    ## if check_option(varargin,'maxDist'):
    if varargin == 'maxDist':                   # Unsure on this translation
        # Potentially ebsd should be a class with how it appears structured here
        ##   xyDist = sqrt((ebsd.prop.x(Dl)-ebsd.prop.x(Dr)).^2 + (ebsd.prop.y(Dl)-ebsd.prop.y(Dr)).^2)
        xyDist = np.sqrt((ebsd.prop.x(Dl)-ebsd.prop.x(Dr)) ** 2 + (ebsd.prop.y(Dl)-ebsd.prop.y(Dr)) ** 2)
        ## dx = sqrt(sum((max(ebsd.unitCell)-min(ebsd.unitCell)).^2))
        dx = np.sqrt(sum((max(ebsd.unitCell) - min(ebsd.unitCell)) ** 2))
        ## maxDist = get_option(varargin,'maxDist',3*dx)                  # Skipping this line for the time being
        ## maxDist = get_option(varargin,'maxDist',inf)                   # This was commented out from mtex
    else:
        maxDist = 0

    ## connect = zeros(size(Dl))
    connect = np.zeros(np.shape(Dl))

    ## for p = 1:numel(ebsd.phaseMap)
    for p in 1:np.size(ebsd.phaseMap)
        # neighboured cells Dl and Dr have the same phase
        if maxDist > 0:
            ## ndx = ebsd.phaseId(Dl) == p & ebsd.phaseId(Dr) == p & xyDist < maxDist       # returns index if all true
            if p == ebsd.phaseId(Dl) and p == ebsd.phaseId(Dr) and xyDist < maxDist:
                ndx = p
        else:
        ## ndx = ebsd.phaseId(Dl) == p & ebsd.phaseId(Dr) == p                          # returns index if all true
            if p == ebsd.phaseId(Dl) and p == ebsd.phaseId(Dr):
                ndx = p

        ## connect(ndx) = true                                                            # Can't find connect function

        # check, whether they are indexed
        # ndx = ndx & ebsd.isIndexed(Dl) & ebsd.isIndexed(Dr)                             # returns index if all true

        # now check for the grain boundary criterion
        if any(ndx):
            ## connect(ndx) = feval(['gbc_' gbc], ebsd.rotations,ebsd.CSList{p},Dl(ndx),Dr(ndx),gbcValue{p},varargin{:})
            connect(ndx) = gbc_angle(ebsd.rotations,ebsd.CSList[p],Dl(ndx),Dr(ndx),gbcValue[p],varargin[:])

    # adjacency of cells that have no common boundary
    ind = connect > 0
    # A_Do = sparse(double(Dl(ind)),double(Dr(ind)),connect(ind),length(ebsd),length(ebsd))
    A_Do = sparse.spmatrix(float(Dl(ind)), float(Dr(ind)), connect(ind), len(ebsd), len(ebsd))
    if check_option(varargin,'mcl'):
        param = get_option(varargin,'mcl')
        if isempty(param), param = 1.4
        if max(param.shape) == 1, param = [param,4]
        A_Do = mclComponents(A_Do,param(1),param(2))
        A_Db = sparse.spmatrix(float(Dl), float(Dr), true, max(ebsd.shape), max(ebsd.shape)) and not A_Do
    else:
        A_Db = sparse.spmatrix(float(Dl(connect<1)), float(Dr(connect<1)), true, max(ebsd.shape), max(ebsd.shape))

    A_Do = A_Do | A_Do.'

    # adjacency of cells that have a common boundary
    A_Db = A_Db | A_Db.'

    # compute I_DG connected components of A_Do
    # I_DG - incidence matrix cells to grains
    ## I_DG = sparse(1:max(ebsd.shape), float(connectedComponents(A_Do)),1)
    I_DG = sparse.spmatrix(1:max(ebsd.shape), float(sparse.csgraph.connected_components(A_Do)),1)

    return A_Db, I_DG

if __name__ == '__main__':
    # Determine which cells to connect
    A_Db, I_DG = do_segmentation(I_FD, ebsd, varargin[:])
    # A_Db - neighboring cells with grain boundary
    # I_DG - incidence matrix cells to grains
