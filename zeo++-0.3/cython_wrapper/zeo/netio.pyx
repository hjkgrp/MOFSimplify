# distutils: language = c++
# distutils: sources = ../networkio.cc

from netstorage cimport AtomNetwork, VoronoiNetwork
from netstorage cimport ATOM_NETWORK, VORONOI_NETWORK
# Define the python definitions for the zeo++ functions

# Easier to implement in python
#def void parseFilename(const char* filename, char* name, char* extension):

#def bint checkInputFile(char* filename)

def readCiffile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not  readCIFFile(c_filename, atmnet.thisptr, radialflag):
        raise ValueError        # Find the appropriate error and return it
    return atmnet

def readArcfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readARCFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readCucfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readCUCFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readCssrfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readCSSRFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readV1file(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readV1File(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def writeCssrfile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToCSSR(c_filename, c_atmnet):
        raise IOError

def writeCiffile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToCIF(c_filename, c_atmnet):
        raise IOError

def writeV1file(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToV1(c_filename, c_atmnet):
        raise IOError

def writeNt2file(filename, vornet, minRad = None):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    if minRad:
        if not writeToNt2(c_filename, c_vornet_ptr, minRad):
            raise IOError
    else:
        if not writeToNt2(c_filename, c_vornet_ptr):
            raise IOError


def writeXyzfile(filename, atmnet, supercell_flag, is_duplicate_perimeter_atoms):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToXYZ(c_filename, c_atmnet, supercell_flag, 
            is_duplicate_perimeter_atoms):
        raise IOError

def writeVtkfile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToVTK(c_filename, c_atmnet):
        raise IOError

def writeMopacfile(filename, atmnet, supercell_flag):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToMOPAC(c_filename, c_atmnet, supercell_flag):
        raise IOError


