# distutils: language = c++
# distutils: sources = ../graphstorage.cc
"""
Cython declarations file for Zeo++ graphstorage module.
Declares Zeo++ DIJSTRA_NODE, DIJKSTRA_NETWORK classes and the associated
python wrappers, DijkstraNode, DijkstraNetwork
"""

__author__ = "Bharat Kumar Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "Jan 6, 2014"

from netstorage cimport VORONOI_NETWORK

cdef extern from "../../graphstorage.h":
    cdef cppclass DIJKSTRA_NODE:
        int id
        double x, y, z
        bint active
        DIJKSTRA_NODE() except +
    cdef cppclass DIJKSTRA_NETWORK:
        DIJKSTRA_NETWORK() except +
cdef extern from "../../graphstorage.h" namespace "DIJKSTRA_NETWORK":
    cdef void buildDijkstraNetwork(VORONOI_NETWORK*, DIJKSTRA_NETWORK*)


#cdef class DijkstraNode:
#    """
#    Cython wrapper class for Zeo++ DIJKSTRA_NODE class.
#    """
#    cdef DIJKSTRA_NODE* thisptr


cdef class DijkstraNetwork:
    """
    Cython wrapper class for Zeo++ DIJKSTRA_NETWORK class.
    """
    cdef DIJKSTRA_NETWORK* thisptr
