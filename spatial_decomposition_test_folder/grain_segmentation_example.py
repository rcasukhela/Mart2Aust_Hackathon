# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 09:22:49 2022

@author: agerlt
"""

# Imports
import tempfile
from diffpy.structure import Atom, Lattice, Structure
import numpy as np
import networkx as nx
from numpy import pi as pi
from scipy import sparse
from orix import data, io, plot
from orix.crystal_map import CrystalMap, Phase, PhaseList
from orix.quaternion import Orientation, Rotation, symmetry
from orix.vector import Vector3d
import matplotlib.pyplot as plt

# ============ #
#   Functions  #
# ============ #


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

# ============ #
#     Code     #
# ============ #

# make a file to download the practice dataset into
tempdir = tempfile.mkdtemp() + "/"
# download a dataset to look at
# Download it, then load it in as an crystal map object
xmap = data.sdss_ferrite_austenite(allow_download=True)
# pick a cutoff angle (ie, disorientation angle)
cutoff_angle = np.deg2rad(15)

# Assuming regular grid, MTEX's grain segmentation is super trivial. First
# generate the adjacency matrix (this next chunk will get replaced by 
# spatial decomp for non-regular grids, but for now, that doesn't matter)

id_map = np.arange(xmap.size).reshape(xmap.shape)
lr_angle = ((xmap[1:, :].rotations)*~(xmap[:-1, :].rotations)).angle
ud_angle = ((xmap[:, 1:].rotations)*~(xmap[:, :-1].rotations)).angle
lr_all_edges = np.vstack([id_map[:-1, :].flatten(), id_map[1:, :].flatten()]).T
ud_all_edges = np.vstack([id_map[:, :-1].flatten(), id_map[:, 1:].flatten()]).T
lr_edges = lr_all_edges[lr_angle < cutoff_angle]
ud_edges = ud_all_edges[ud_angle < cutoff_angle]
all_internal_edges = np.vstack([lr_edges, ud_edges])
all_external_edges = np.vstack([lr_all_edges[lr_angle > cutoff_angle],
                                lr_all_edges[lr_angle > cutoff_angle]])
A_internal = (1+all_internal_edges[:, 0]*0,
              (all_internal_edges[:, 0], all_internal_edges[:, 1]))
A = sparse.coo_matrix(A_internal, shape=(xmap.size, xmap.size))
# the A above is the adjacency where connections that are above the cutoff are
# gone. Stick that into networkx, because they already have the "grouping"
# algorithm and its super fast
G = nx.from_scipy_sparse_array(A)
feature_clusters = list(nx.connected_components(G))
grain_ids = np.zeros(xmap.size)
for i in range(len(feature_clusters)):
    grain_ids[list(feature_clusters[i])] = i
grain_ids_square = grain_ids.reshape(xmap.shape)
# feature_clusters is all you need to make shapely objects. it is a list of
# all the pixel ids that go to every grain. grain_ids is the grain number for
# every pixel, and grain_ids_square is just that as a 2D array, so you can
# look at it more easily

# Visualization stuff:
plt.close('all')
plt.figure()
xmap.plot()
plt.figure()
plt.imshow(grain_ids_square % pi)

# FROM HERE: you SHOULD be able to make shapely objects from each of the lists
# stored in feature_clusters (its a list of lists on integers i, where V[i] 
# gives the coordinates of the two ends of the associated edge)

# ask me if you need clarification, this is very hand-wavey
