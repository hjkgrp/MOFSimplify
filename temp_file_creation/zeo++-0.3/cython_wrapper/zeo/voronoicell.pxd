"""
Cython declarations file for Zeo++ classes in geometry module.
Declares Zeo++ XYZ and Point Classes and the associated Python wrappers.
"""
__author__ = "Bharat Medasani"
__date__ = "2014-02-17"

from libcpp.map cimport map as cmap
from libcpp.set cimport set as cset
from libcpp.vector cimport vector

from zeo.geometry cimport CPoint
from zeo.netstorage cimport ATOM_NETWORK, VORONOI_NETWORK

cdef extern from "../../voronoicell.h":
    cdef cppclass VOR_FACE:
        VOR_FACE(vector[CPoint], ATOM_NETWORK*, VORONOI_NETWORK*) except +
        vector[CPoint] vertices "orderedVertices"
        vector[int] node_ids "nodeIDs"


    cdef cppclass VOR_CELL:
        VOR_CELL() except +
        vector[VOR_FACE] faces
        int num_vertices "numVertices"
        #cmap[CPoint, int, bint
        cmap[int,int] id_mappings "idMappings"
        cmap[int, vector[int]] reverse_id_mappings "reverseIDMappings"
        cmap[int, CPoint] vertex_coords "vertexCoords"
        vector[cset[int]] edge_connections "edgeConnections"

        void add_edge "addEdge" (CPoint, CPoint to)
        void add_face "addFace" (VOR_FACE)

    cdef cppclass BASIC_VCELL:
        BASIC_VCELL() except +


cdef class VorFace:
    cdef  VOR_FACE* thisptr

cdef class VorCell:
    cdef VOR_CELL* thisptr

cdef class BasicVCell:
    cdef BASIC_VCELL* thisptr
