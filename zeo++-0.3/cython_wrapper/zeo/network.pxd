# distutils: language = c++
# distutils: sources = ../network.cc

from libcpp.vector cimport vector 
from libcpp.string cimport string
from zeo.netstorage cimport ATOM_NETWORK, VORONOI_NETWORK
from zeo.voronoicell cimport VOR_CELL, BASIC_VCELL


# At present  the return value of performVoronoiDecomp is void*
# Compile it after void* is changed to bool in the original source file
cdef extern from "../../network.h":
    cdef bint performVoronoiDecomp(bint, ATOM_NETWORK*, VORONOI_NETWORK*, 
            vector[VOR_CELL]*, bint, vector[BASIC_VCELL]*)

    cdef void calculateFreeSphereParameters(VORONOI_NETWORK*, char*, bint)

    cdef void viewVoronoiDecomp(ATOM_NETWORK*, double, string)

    cdef void loadRadii(ATOM_NETWORK*)

    cdef void loadMass(bool, ATOM_NETWORK*)

cdef extern from "../../area_and_volume.h":
    cdef void visVoro(char* name, double probeRad, int skel_a, int skel_b, int skel_c,
            VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet)
