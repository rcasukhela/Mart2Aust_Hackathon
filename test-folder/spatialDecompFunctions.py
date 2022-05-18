#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:00:47 2022

@author: paytone
"""

def gbc_angle(q, CS, D_l, D_r, threshold=5.):
    '''
    THIS IS AN INITIAL DRAFT TRANSLATION OF GMTEX's BC_ANGLE
    
    
    # the inputs are: quaternion(ebsd.rotations),ebsd.CSList{p},Dl(ndx),Dr(ndx),gbcValue(p),varargin{:})

    
    #% FOR TESTING
    ################## load materials and image 
    # this whole section will be replaced by graph cut function eventually
    path = r'/Users/paytone/Documents/GitHub/maxflow_for_matsci/Data/steel_ebsd.ang'
    
    # Read each column from the file
    euler1, euler2, euler3, x, y, iq, dp, phase_id, sem, fit  = np.loadtxt(path, unpack=True)
    
    phase_id = np.int32(phase_id)
    
    # Create a Rotation object from Euler angles
    euler_angles = np.column_stack((euler1, euler2, euler3))
    rotations = Rotation.from_euler(euler_angles)
    
    # Create a property dictionary
    properties = dict(iq=iq, dp=dp)
    
    # Create unit cells of the phases
    structures = [
        Structure(
            title="ferrite",
            atoms=[Atom("fe", [0] * 3)],
            lattice=Lattice(0.287, 0.287, 0.287, 90, 90, 90)
        ),
    ]
    phase_list = PhaseList(
        names=["ferrite"],
        point_groups=["432"],
        structures=structures,
    )
    
    # Create a CrystalMap instance
    xmap2 = CrystalMap(
        rotations=rotations,
        phase_id=phase_id,
        x=x,
        y=y,
        phase_list=phase_list,
        prop=properties,
    )
    xmap2.scan_unit = "um"
    
    import numpy
    
    #% get pairs of neighbouring cells {D_l,D_r} in A_D
    #A_D = I_FD'*I_FD==1;
    #[Dl,Dr] = find(triu(A_D,1));
    '''
        
    import numpy as np
    from orix.quaternion import Orientation, Rotation, symmetry, Misorientation
    
    # convert threshold into radians
    threshold = threshold * np.pi / 180.
    
    
    # now check whether the have a misorientation heigher or lower than a threshold
    o1 = q[D_l] # orientation of node u
    o2 = q[D_r] # orientation of node v
    m = Misorientation(o1 * o2.conj) # misorientations between every u and v
    m.symmetry = (CS, CS) # Oh is symmetry (need to un-hard code)
    m = m.map_into_symmetry_reduced_zone() #TODO: The result of this doesn't actually make sense. Why are there two rows?
    
    criterion = np.abs(m[0,:].angle) > np.cos(threshold/2.0)
    
    # Part of mtex code; purpose not clear:
    #if np.any(~criterion):
    #  qcs = symm_l.proper_subgroup #TODO: Look into whether this actually makes sense, would be best to use the same symmetry for both L and R but it wasn't clear how to do this in orix at this time
    #  criterion[~criterion] = np.amax(np.abs(Rotation.dot_outer(m[~criterion],qcs)), axis=1) > np.cos(threshold/2.);
    
    # Here is how mtex gets criterion using orientations, not yet possible in orix
    # o_Dl = orientation(q(Dl),CS,symmetry);
    # o_Dr = orientation(q(Dr),CS,symmetry);
    # criterion = dot(o_Dl,o_Dr) > cos(threshold/2);
    
    return criterion


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
    # H = @(v) eye(3) - 2./(v(:)'*v(:))*(v(:)*v(:)') ;
    v = np.atleast_2d(v)
    H = np.eye(3) - 2. / np.matmul(v, v.T) * np.matmul(v.T, v)
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
    from orix.utilities.utilities import sortrows
    from orix.utilities import utilities

    xy = X
    unitCell = unit_cell

    def clockwise_sort(points_list):
        import math

        def angle_to(point):
        # assumes already mapped to (0,0) center
        # modify as needed to get result in desired range
            return math.atan2(point[1],point[0])

        sorted_points = sorted(points_list, key=angle_to, reverse=True)
        sorted_points = np.roll(sorted_points, 1, axis=0)

        return sorted_points

    unitCell = clockwise_sort(unitCell)

    # add unitCell values to data.
    def add_unitcell_values(X):
        X_new = []
        for elem1 in unitCell:
            for elem2 in X:
                X_new.append(list(elem1+elem2))

        X_new = np.array(X_new)
        return X_new
        
    X_new = add_unitcell_values(X)

    x = X_new[:, 0]
    y = X_new[:, 1]
    #X_new = None

    # remove equal points
    eps = np.amin([np.sqrt(np.diff(unitCell[:,0])**2 + np.diff(unitCell[:,1])**2)])/10.;

    verts = np.squeeze(np.round([x - np.amin(x), y - np.amin(y)]/eps)).T

    v, m, n = utilities.uniquerows(verts)
    v, order = sortrows(v, return_order=True)

    x = np.squeeze(x)
    y = np.squeeze(y)
    m = np.squeeze(m)


    v = np.vstack([x[m][order], y[m][order]]).T

    # set faces
    faces = np.reshape(order[n], [-1,np.size(unitCell[:,1])])
    
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