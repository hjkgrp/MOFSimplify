# distutils: language = c++
# distutils: sources = ../networkstorage.cc
"""
Cython decloarations file for Zeo++ network storage section.
Declares Zeo++ ATOM_NETWORK, VORONOI_NETWORK classes and the associated 
python wrappers AtomNetwork, VoronoiNetwork.
"""
__author__ = "Bharat Kumar Medasani"
__date__ = "2013-12-09"

from libcpp.vector cimport vector

from zeo.geometry cimport CPoint
from zeo.voronoicell cimport VOR_CELL, BASIC_VCELL


cdef extern from "../../networkstorage.h":
    cdef cppclass ATOM:
        ATOM() except +
        double x, y, z
        double radius
        #string type
        #int specialID
        double mass
        double charge 

    cdef cppclass ATOM_NETWORK:
        ATOM_NETWORK() except +
        void copy(ATOM_NETWORK*)
        double a, b, c      # Lattice parameters
        double alpha, beta, gamma   # lattice angles
        int no_atoms "numAtoms"
        vector[ATOM] atoms
        CPoint abc_to_xyz(double, double, double)
        CPoint abc_to_xyz(CPoint)
        CPoint xyz_to_abc(double, double, double)
        CPoint xyz_to_abc(CPoint)

    cdef cppclass VOR_NODE:
        VOR_NODE() except +
        VOR_NODE(double, double, double, double, vector[int])
        double x, y, z
        double rad_stat_sphere

    cdef cppclass VOR_EDGE:
        VOR_EDGE() except +
        VOR_EDGE(int, int, double, int, int, int, double)
        int origin "from" 
        int ending "to"
        double rad_moving_sphere
        double length
        int delta_uc_x, delta_uc_y, delta_uc_z


    cdef cppclass VORONOI_NETWORK:
        VORONOI_NETWORK() except +
        VORONOI_NETWORK prune(double)
        vector[VOR_NODE] nodes
        vector[VOR_EDGE] edges

    cdef bint c_substituteAtoms "substituteAtoms"(ATOM_NETWORK*, ATOM_NETWORK*,
            bint, int*, bint)

    cdef bint c_fracSubstituteAtoms "fracSubstituteAtoms"(ATOM_NETWORK*, 
            ATOM_NETWORK*, bint, double, 
            int, int*, double*, bint)

cdef extern from '../../networkio.h':
    cdef void parseFilename(const char* filename, char* name, char* extension)

    cdef bint checkInputFile(char* filename)

    cdef bint readCIFFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readARCFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readCUCFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readCSSRFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readV1File(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint writeToCSSR(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToCIF(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToV1(char * filename, ATOM_NETWORK *cell)

    cdef bint writeToNt2(char *filename, VORONOI_NETWORK *vornet, double minRad)

    cdef bint writeToNt2(char *filename, VORONOI_NETWORK *vornet)

    cdef bint writeToXYZ(char *filename, ATOM_NETWORK *cell, bint is_supercell,
                         bint is_duplicate_perimeter_atoms)

    cdef bint writeToVTK(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToMOPAC(char *filename, ATOM_NETWORK *cell, bint is_supercell)

    cdef bint writeVornetToXYZ "writeToXYZ"(char *filename, VORONOI_NETWORK*, double)

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

cdef class Atom:
    """
    Cython wrapper class for Zeo++ ATOM class.
    """
    cdef ATOM* thisptr

cdef class AtomNetwork:
    """ 
    Cython wrapper class for Zeo++ ATOM_NETWORK class.
    Contains a pointer to ATOM_NETWORK and a flag denoting whether radius
    for each atomic species is non-zero.
    """
    cdef ATOM_NETWORK* thisptr
    cdef bint rad_flag

cdef class VoronoiNode:
    """
    Cython wrapper class for Zeo++ VOR_NODE class.
    """
    cdef VOR_NODE* thisptr

cdef class VoronoiNetwork:
    """ 
    Cython wrapper class for Zeo++ VORONOI_NETWORK class.
    """
    cdef VORONOI_NETWORK* thisptr
