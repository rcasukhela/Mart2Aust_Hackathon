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

    
    #%% FOR TESTING
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
    m = Misorientation([o1.data, o2.data]) # misorientations between every u and v
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
