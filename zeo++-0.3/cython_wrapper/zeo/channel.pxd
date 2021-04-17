from libcpp.vector cimport vector
from netstorage cimport VORONOI_NETWORK
from graphstorage cimport DIJKSTRA_NETWORK


cdef extern from "../../channel.h":
    cdef cppclass CHANNEL:
        CHANNEL() except +

cdef extern from "../../channel.h" namespace "CHANNEL":
    cdef c_findChannelsInDijkstraNet "findChannels"(DIJKSTRA_NETWORK*, 
            vector[bint] *, vector[CHANNEL] *)
    cdef c_findChannelsInVorNet "findChannels"(VORONOI_NETWORK*, double, 
            vector[bint] *, vector[CHANNEL] *)

cdef class Channel:
    cdef CHANNEL* thisptr
