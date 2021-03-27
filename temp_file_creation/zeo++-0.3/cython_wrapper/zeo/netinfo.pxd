# distutils: language = c++
# distutils: sources = ../networkinfo.cc

from libcpp.map cimport map
from libcpp.string cimport string

cdef extern from "../../networkinfo.h":
    cdef void zeo_initializeRadTable "initializeRadTable"()

    cdef void zeo_initializeCovRadTable "initializeCovRadTable"()

    cdef void zeo_initializeMassTable "initializeMassTable"()

    cdef void zeo_initializeAtomCharacterTable "initializeAtomCharacterTable"()

    cdef void zeo_initializeAtomicNumberTable "initializeAtomicNumberTable"()

    cdef void zeo_readRadTable "readRadTable"(char *filename)

    cdef void zeo_readMassTable "readMassTable"(char *filename)

    cdef double zeo_lookupRadius "lookupRadius"(string atomType, bint radial)

    cdef double zeo_lookupCovRadius "lookupCovRadius"(string atomType)

    cdef double zeo_lookupMass "lookupMass"(string atomType)

    cdef int zeo_lookupAtomicNumber "lookupAtomicNumber"(string atomType)

    cdef bint zeo_isMetal "isMetal"(string atomType)

