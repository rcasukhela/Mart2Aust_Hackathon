#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:00:47 2022

@author: paytone
"""
import numpy as np
from scipy import sparse



def basic_spatial_decomp(xmap, unit_cell=None):
    """
    converts a crystal map into a series of cells (D) with vertices (V),
    where every cell correlates with a single pixel, and forms a convex hull
    surrounding the voronoi area around each voxel in the crystal map.

    Parameters
    ----------
    xmap: orix.xmap object
        the 2D or 3D ebsd map of orientation and phase data

    unit_cell: None or Voronoi
        If 'None', will calculate the expected cells and vertexes assuming a
        repeating rectangular grid. (hex is currently unsupported)
        If 'Voronoi', will calculate V and D using qhull's voronoi tesselation

    Returns
    -------
    V:
        n x m numpy array of vertex coordinates
    D:
        4 x m numpy array of vertex indicies that form D. the coordinates of
        the vertexes of D can be found as V[D,:], or the Vertexes of the xth
        entry in D as V[D[x],:]
    """
    if unit_cell is not None:
        raise NotImplemented("bingus bongus")
    xx, yy = np.meshgrid(
        np.arange(xmap.shape[1] + 1)*xmap.dx - xmap.dx/2,
        np.arange(xmap.shape[0] + 1)*xmap.dy - xmap.dy/2
        )
    V = np.vstack([xx.flatten(), yy.flatten()]).T
    N = np.arange(xx.size).reshape(xx.shape)
    D = np.vstack([
        N[:-1, :-1].flatten(),
        N[:-1, 1:].flatten(),
        N[1:, :-1].flatten(),
        N[1:, 1:].flatten()
        ]).T
    big_F = np.vstack([D[:, (0, 1)], D[:, (0, 2)], D[:, (1, 3)], D[:, (2, 3)]])
    F, ie = np.unique(big_F, axis=0, return_inverse=True)
    F_ID = sparse.coo_matrix(
        (ie*0 + 1, (np.arange(xmap.size).repeat(4), ie)),
        shape=(xmap.size, F.size)
        )
    return(V, D, F_ID)


def erase_linearly_dependent_points(points):
    '''
    subfunction to remove linearly dependent points.

    Inputs:
    --------------
    k : ???
        ???

    Outputs:
    --------------
    ??? : ???
        ???

    Dependencies:
    --------------
    from scipy.spatial import ConvexHull
    '''
    import numpy as np
    from scipy.spatial import ConvexHull
    k = ConvexHull(points)

    # erase all linear dependent points
    angle = np.arctan2( 
        x[k.vertices[0:-1]]-x[k.vertices[1:]],
        y[k.vertices[0:-1]]-y[k.vertices[1:]]
    )
    test = np.abs(np.diff(angle))>np.spacing(1.0)
    k2 = k.vertices[np.concatenate([[True], test, [True]])]
    boundingX = [x[k2], y[k2]]

    cellRot=0

################################################################################
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
                    # REQUIRES SUBROUTINE CONSTRUCTION #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
    
    unitCell = util.regularPoly(4,xstp,cellRot)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
                    # REQUIRES SUBROUTINE CONSTRUCTION #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
################################################################################

    print(unitCell)
    radius = np.mean( np.sqrt(np.sum(unitCell**2,2)) )
    edgeLength = np.sqrt( np.sum( np.diff(boundingX)**2, axis = 2) )

    return k