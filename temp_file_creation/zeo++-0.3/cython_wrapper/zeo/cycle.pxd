"""
Cython declarations file for Zeo++ cycle.h section
Declares Zeo++ CYCLE class and the associated python wrapper, Cycle.
"""

__author__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "Jan 6, 2014"

from libcpp.vector cimport vector

from zeo.netstorage cimport VORONOI_NETWORK, ATOM_NETWORK
from zeo.geometry cimport XYZ
from zeo.graphstorage cimport DIJKSTRA_NODE

cdef extern from "../../cycle.h":
    cdef cppclass CYCLE:
        CYCLE() except +
        double length
        vector[DIJKSTRA_NODE] nodes

    cdef bint compute_4cycle(VORONOI_NETWORK*, vector[CYCLE]*, bint, int)

    cdef void centroid(CYCLE*, XYZ*, vector[int]*)

    cdef void face_center(ATOM_NETWORK*, vector[XYZ]*)


cdef class Cycle:
    """
    Cython wrapper class for Zeo++ CYCLE class.
    Contains a pointer to CYCLE
    """
    cdef CYCLE* thisptr

