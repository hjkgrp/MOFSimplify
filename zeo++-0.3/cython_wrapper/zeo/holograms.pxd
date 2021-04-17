from libcpp.vector cimport vector
from netstorage cimport VORONOI_NETWORK

cdef extern from "../../holograms.h":
    cdef void analyze_accessible_voronoi_pre_segment(VORONOI_NETWORK *vornet, 
            float probeRad, vector[bint] *accessInfo, char *name, 
            char *bin_directory) 
