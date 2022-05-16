#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:00:47 2022

@author: paytone
"""

def householderMatrix(v):
    '''
    Parameters
    ----------
    v : TYPE
        DESCRIPTION.

    Returns
    -------
    H : TYPE
        DESCRIPTION.
        
    For Testing:
    H_out = np.loadtxt('householder_output.txt', delimiter=',')
    H_in = np.loadtxt('householder_input.txt', delimiter=',')
    H_test = householderMatrix(H_in)

    '''
    import numpy as np
    v = np.atleast_2d(v).T
    H = np.eye(3) - 2. / np.matmul(v.T, v) * np.matmul(v, v.T)
    return H

def translationMatrix(s):
    '''
    Parameters
    ----------
    v\s : TYPE
        DESCRIPTION.

    Returns
    -------
    T : TYPE
        DESCRIPTION.
        
    For Testing:
    T_out = np.loadtxt('translation_output.txt', delimiter=',')
    T_in = np.loadtxt('translation_input.txt', delimiter=',')
    T_test = translationMatrix(T_in)
    '''
    import numpy as np
    T  = np.array([[ 1., 0., s[0]],[0., 1., s[1]],[0., 0., 1.]])
    return T




def generateUnitCells(xy, unitCell):
    '''
    % generates a list of patches according to spatial coordinates and the unitCell
    %
    % Inputs
    %  xy       - midpoints of the cells
    %  unitCell - spatial coordinates of the unit cell
    %
    % Parameters
    %  v     - list of vertices
    %  faces - list of faces
    
    # For Testing:
    xy = np.loadtxt('spatialDecomposition_input_X.csv', delimiter=',')
    unitCell = np.loadtxt('spatialDecomposition_input_unitCell.csv', delimiter=',') 
    v, faces = generateUnitCells(xy, unitCell)  
    randomface = np.random.randint(np.size(xy[:,0]))
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(v[faces[randomface,:],0], v[faces[randomface,:],1], '--')
    plt.axis('equal')
    plt.show
    '''
    import numpy as np
    
    # compute the vertices
    x = np.reshape(np.tile(np.atleast_2d(xy[:,0]).T, [1, np.shape(unitCell[:,0])[0]]) + np.tile(unitCell[:,0],[np.shape(xy[:,0])[0],1]), [1,-1])
    y = np.reshape(np.tile(np.atleast_2d(xy[:,1]).T, [1, np.shape(unitCell[:,1])[0]]) + np.tile(unitCell[:,1],[np.shape(xy[:,1])[0],1]), [1,-1])
    
    # remove equal points
    eps = np.amin([np.sqrt(np.diff(unitCell[:,0])**2 + np.diff(unitCell[:,1])**2)])/10.;
    
    verts = np.squeeze(np.round([x - np.amin(x), y - np.amin(y)]/eps)).T
    
    from orix.utilities import utilities
    v, m, n = utilities.uniquerows(verts);
    
    x = np.squeeze(x)
    y = np.squeeze(y)
    m = np.squeeze(m)
    
    v = np.vstack([x[m], y[m]]).T
    
    # set faces
    faces = np.reshape(n, [-1,np.size(unitCell[:,1])]);
    
    return v, faces

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
