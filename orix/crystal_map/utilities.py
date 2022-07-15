#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 09:08:49 2022

@author: paytone
"""
import numpy as np

def spatial_decomposition(X, *args):
    '''
    % decomposite the spatial domain into cells D with vertices V,
    %
    % Output
    %  V - list of vertices
    %  F - list of faces
    %  I_FD - incidence matrix between faces to cells

    % compute voronoi decomposition
    % V - list of vertices of the Voronoi cells
    % D   - cell array of Vornoi cells with centers X_D ordered accordingly
    '''
    # Imports
    import getopt
    from numpy import cumsum, zeros, unique, sort
    from spatialDecompFunctions import generateUnitCells
    from scipy.spatial import Voronoi
    from scipy.sparse import csr_matrix, coo_array
    from orix.utilities.utilities import uniquerows

    if unit_cell.all() == None:
        unit_cell = calcUnitCell(X)

    if args[0] == 'unit_cell':

        # compute the vertices
        [V, faces] = generateUnitCells(X, unit_cell, args[:])
        # NOTE: V and faces do not give exact values as compared to the MATLAB
        # implementation. However, Eric and Rohan have confirmed that V and faces
        # at least have the same shape as the MATLAB implementation.

        D = np.empty(len(X),dtype=object)

        for k in range(X.shape[0]):
            D[k] = faces[k, :]

    else:    
        var_arg_in = args[0]
        dummy_coordinates = calcBoundary(X, unit_cell, var_arg_in)

        vor = Voronoi(np.vstack([X, dummy_coordinates]), 
                        qhull_options = 'Q5 Q6 Qs') #,'QbB'

        V = vor.vertices
        D = vor.regions
        D = D[0:X.shape[0]]

    # now we need some adjacencies and incidences
    iv = np.hstack(D)        # nodes incident to cells D
    iid = zeros(len(iv), dtype=np.int64)   # number the cells

    # Some MATLAB stuff goin on here... : p = [0; cumsum(cellfun('prodofsize',D))];
    D_prod = matlab_prod_of_size(D)
    p = np.hstack([0, np.cumsum(D_prod)-1])

    # original matlab: for k=1:numel(D), id(p(k)+1:p(k+1)) = k; end
    for k in range(len(D)):
        iid[p[k]+1 : p[k+1]] = k

    # next vertex
    # indx = 2:numel(iv)+1;
    ind_x = np.arange(1,len(iv)+1)
    # indx(p(2:end)) = p(1:end-1)+1;
    ind_x[p[1:]] = p[0:-1] + 1
    # ivn = iv(indx);
    iv_n = iv[ind_x]

    # edges list
    F = np.vstack([iv[:], iv_n[:]]).T

    # should be unique (i.e one edge is incident to two cells D)
    F, _, ie = uniquerows(F)

    # faces incident to cells, F x D
    #original matlab: I_FD = sparse(ie,id,1); 
    # Matlab stores as csr matrix, so we use this class below
    data = np.ones(np.shape(ie))
    I_FD = csr_matrix( (data, (ie, iid)))
    #I_FD = csr_matrix( (data, (ie, iid)), shape=[ie.shape[0], ie.shape[0]]) # could also use csc_matrix() if it improves
                                # performance !
    '''
    NOTE: may need to use different sparse matrix. The matrix is incorrect as
          compared to the MATLAB implementation.
          
          Could also be causing problems in the Householder matrix initialization.
    '''
    

    return V, F, I_FD

def calcBoundary(X, unit_cell, var_arg_in='hull'):
    '''
    dummy coordinates so that the voronoi-cells of X are finite

    Inputs:
    --------------
    X : n x 2 numpy array of [x,y] vertex coordinates

    unit_cell : n x 2 numpy array of [x,y] coordinates of "pixel" boundary (e.g. hexagon or square)

    var_arg_in : ???
    ???


    Outputs:
    --------------
    ??? : ???
    ???

    Dependencies:
    --------------
    import getopt
    import numpy as np
    from scipy.spatial import ConvexHull
    from numpy import arctan2, diff, cumsum, matrix, squeeze, unique
    from statistics import mean
    from math import sqrt, floor, copysign
    '''
    # Imports
    import getopt
    import numpy as np
    from scipy.spatial import ConvexHull
    from numpy import arctan2, cumsum, matrix, squeeze, unique, linspace
    from statistics import mean
    from math import sqrt, floor, copysign
    from orix.quaternion.rotation import Rotation
    import orix.vector as vect
    from orix.utilities.utilities import uniquerows
    
    from spatialDecompFunctions import householderMatrix, translationMatrix#, erase_linearly_dependent_points

    dummy_coordinates = []

    boundary = 'hull'
    boundary = str(var_arg_in)

    if boundary.isalpha():
        
        if boundary.lower() == 'tight':
            
            k = alphashape.alphashape(X, 0.95) #TODO: 0.95 matches mtex, but this would be a good user variable
            k2 = list(zip(*k.exterior.coords.xy))
            locations = []
            for item in k2:
                locations.append(np.where(np.all(item==X, axis=1)))
            locations = np.squeeze(locations)
            
            bounding_X  = erase_linearly_dependent_points(X, locations)
            
            #Testing
            #import matplotlib.pyplot as plt
            #plt.figure()
            #lt.plot(bounding_X[:,0], bounding_X[:,1])
            #plt.show()

        elif (boundary.lower() == 'hull' or
            boundary.lower() == 'convexhull'):
            
            bounding_X  = erase_linearly_dependent_points(X)
            
            #Testing
            #import matplotlib.pyplot as plt
            #plt.figure()
            #plt.plot(bounding_X[:,0], bounding_X[:,1])
            #plt.show()
            
        elif boundary.lower() == 'cube':
            # set up a rectangular box
            envelope_X = [np.amin(X), np.amax(X)]
            bounding_X = [
                envelope_X[1], envelope_X[3],
                envelope_X[2], envelope_X[3],
                envelope_X[2], envelope_X[4],
                envelope_X[1], envelope_X[4],
                envelope_X[1], envelope_X[3],
                ]

        else:
            raise ValueError('Unknown boundary type. Available options are \
            ''hull'', ''convexhull'' and ''cube''.')

    elif isinstance(boundary, float):
        bounding_X = boundary
    
    radius = np.mean(np.sqrt(np.sum( (unit_cell**2), axis=1)))
    
    edge_length = ((matlabdiff(bounding_X)**2).sum(axis=1, keepdims=True))**0.5 
    #edgeLength = sqrt(sum(diff(boundingX).^2,2));
    
    # fill each line segment with nodes every 20 points (in average)
    nto = np.int32(np.fix( (edge_length>0)*4 ))

    #nto = floor( edge_length*(2*radius) )
    #nto = fix((edgeLength>0)*4); fix(edgeLength*(2*radius));

    csum = np.int32(cumsum(np.vstack([0, nto+1])))

    bounding_X = left_hand_assignment(bounding_X, csum)
    
    # interpolation
    for k in range(nto.shape[0]):
        for dim in range(nto.shape[1]+1):
            bounding_X[csum[k]:csum[k+1]+1, dim] = linspace(
                bounding_X[csum[k],dim],
                bounding_X[csum[k+1], dim],
                nto[k,0]+2
            )

    # homogeneous coordinates
    X = np.hstack([X, np.ones([X.shape[0],1])])
    bounding_X = np.hstack([bounding_X, np.ones([bounding_X.shape[0],1])])

    # direction of the edge
    edge_direction = matlabdiff(bounding_X)
    edge_angle = arctan2( edge_direction[:, 1], edge_direction[:,0] )
    edge_length = np.sqrt( (edge_direction**2).sum(axis = 1, keepdims=True))

    # shift the starting vertex
    r = Rotation.from_axes_angles(vect.Vector3d([0, 0, 1]), edge_angle)
    a = vect.Vector3d([0, radius, 1])
    b_X = np.squeeze(r * a)
    offset_X = np.squeeze(b_X.data) - np.squeeze(bounding_X[0:-1, :])
 
    dummy_coordinates = np.array([[0,0]])
    for k in range(bounding_X.shape[0]-1):

        # mirror the point set X on each edge
        p_X = np.matmul(X,
            -1*( 
                np.matmul(
                    np.matmul(
                        translationMatrix(offset_X[k, :]),
                        householderMatrix(edge_direction[k, :])
                    ), 
                    translationMatrix(offset_X[k, :])
                ).T
            )
        )

        '''
        NOTE FROM ERIC AND ROHAN: We are having problems with getting the proper
        dummy_coordinates. The bounding box only has 2 of 4 lines, we believe
        this is where the problem is. MATLAB is doing bizarre things with the
        element-wise division!!!!! We need to contact somebody about this behavior.
        Python is doing what we expect but it's not the right answer.
        
        UPDATE FROM ERIC: Changing the order of operations in Python gives the same result as Matlab.
        '''

        # distance between original and mirrored point
        dist = np.sqrt(np.sum((X[:,0:-1]-p_X[:,0:-1])**2, axis=1))

        intend_X = 2. * radius * np.sign(edge_direction[k,0:2])

        # now try to delete unnecessary points
        m = 2        
        
        # NEED TO REBUILD LOOP STRUCTURE HERE FOR PARITY - JC
        #tmp_X = p_X[np.argwhere(dist < m*radius), 0:-1]
        
        while(True):
        
            tmp_X = p_X[np.argwhere(dist < m*radius), 0:2]
            
            #right = (bsxfun(@minus, tmpX, boundingX(k,1:2)  - intendX ) * edgeDirection(k,1:2)') < 0;
            #left  = (bsxfun(@minus, tmpX, boundingX(k+1,1:2)+ intendX ) * edgeDirection(k,1:2)') > 0;
   
            right = np.matmul(tmp_X - np.tile( bounding_X[k, 0:2]   - intend_X, [np.shape(tmp_X)[0], 1]), edge_direction[k, 0:2].T) < 0
            left  = np.matmul(tmp_X - np.tile( bounding_X[k+1, 0:2] + intend_X, [np.shape(tmp_X)[0], 1]), edge_direction[k, 0:2].T) > 0

            tmp_X = tmp_X[~np.any(np.vstack([right, left]), axis=0), 0:2]

            if edge_length[k] / tmp_X.shape[0] < radius/3.:
                break
            elif m < 2.**7:
                m = m * 2
            elif m < 2.**7 + 100.0:
                m = m + 10.0
            else:
                break
        
        print(np.shape(tmp_X))
        print(np.shape(dummy_coordinates))

        if tmp_X.size > 0:
            dummy_coordinates = np.vstack([dummy_coordinates, tmp_X])

    dummy_coordinates, _, _ = uniquerows(dummy_coordinates[1:,:])
    
    return dummy_coordinates

def matlab_prod_of_size(D):
    '''
    Testing:
    test1 = [[1, 2, 3], [1, 2], [1]]
    matlab_prod_of_size(test1)
    '''
    prod_of_size = []
    for elem in D:
        prod_of_size.append(len(elem))
    return prod_of_size

def erase_linearly_dependent_points(X, k):
    '''
    subfunction to remove linearly dependent points.

    Inputs:
    --------------
    X : n x 2 numpy array of coordinates
    k : list of coordinate pairs corresponding to the convex or concave hull

    Outputs:
    --------------
    boundingX : n x 2 numpy array for bounding box

    Dependencies:
    --------------
    from scipy.spatial import ConvexHull
    '''
    from orix.utilities.utilities import regularPoly
    import numpy as np
    from scipy.spatial import ConvexHull
    from orix.utilities.utilities import sortrows
    import alphashape

    # erase all linear dependent points
    angle = np.arctan2(X[k[0:-1],0]-X[k[1:],0],
        X[k[0:-1],1]-X[k[1:],1])
    test = np.abs(np.diff(angle))>np.spacing(1.0)
    k2 = k[np.concatenate([[True], test, [True]])]
    
    # Build the coord set of the convex hull
    # Pad on the end->start connection to close hull
    # If the origin is contained as a boundary point
    # then we need to shift the array values so the origin
    # is the end->start wrap
    boundingX = X[k2,:]
    
    zero_row = -1
    for i in range(boundingX.shape[0]):
        if np.all(boundingX[i]==0):
            zero_row = i
    
    if (zero_row != -1):
        boundingX = np.roll(boundingX, -zero_row, axis=0)
    
    boundingX = np.vstack([boundingX , boundingX[0,:] ])
 
    return boundingX

def left_hand_assignment(X, a):
    '''
    Attempts to replicate MATLAB left-hand assignment for a 2D array.
 

    Inputs:
    --------------
    X : 2D numpy array
 

    Outputs:
    --------------
    X : 2D numpy array
        This is different from the X passed in, it will return a larger array with
        more zeros.
 
    Dependencies:
    --------------
    import numpy as np
    '''
 
    import numpy as np
    import warnings
    
    if a.dtype != 'int':
        warnings.warn('parameter ''a'' must be of integer type. Converting ''a'' into integers and moving on...')

    a = np.int32(a)

    bound_X_init = X
    ncols = X.shape[1]
 
    max_a = np.amax(a)
    bound_X_fin = np.zeros([max_a+1, ncols])
  
    bound_X_fin[0:bound_X_init.shape[0], :] = bound_X_init

    for i, elem in enumerate(a[0:-1]):
        bound_X_fin[elem, :] = bound_X_init[i, :]

    return bound_X_fin

def matlabdiff(myArray):
    return myArray[1:,:] - myArray[0:-1,:]

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

