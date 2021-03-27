"""
Cython declarations for Zeo++ area_and_volume calculations.
Declares functions for calculating the surface area and volume
of channels and pockets.
"""
from libcpp.string cimport string
from zeo.netstorage cimport ATOM_NETWORK

cdef extern from "../../area_and_volume.h":
    cdef string calcAV(ATOM_NETWORK*, ATOM_NETWORK*, bint, double, double, 
            int, bint, double, double)

    cdef string calcASA(ATOM_NETWORK*, ATOM_NETWORK*, bint, double, double,
            int, bint, bint)

    cdef double calcDensity(ATOM_NETWORK*) 

