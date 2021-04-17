from libcpp.vector cimport vector

#from zeo.geometry cimport CPoint
#from zeo.netstorage cimport AtomNetwork, VoronoiNetwork
#from zeo.netstorage cimport ATOM_NETWORK, VORONOI_NETWORK

#cdef class VorFace:
#    #cdef VOR_FACE *thiptr
#    def __cinit__(self, vertices,  atmnet, vornet):
#        cdef vector[CPoint] c_vertices = (<vertices
#        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
#        cdef VORONOI_NETWORK* c_vornetptr = (<VoronoiNetwork?>vornet).thisptr
#        self.thisptr = new VOR_FACE(c_vertices, c_atmnetptr, c_vornetptr)
#
#    def __dealloc__(self):
#        del self.thisptr

cdef class VorCell:
    #cdef VOR_CELL *thiptr
    def __cinit__(self):
        self.thisptr = new VOR_CELL()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

cdef class BasicVCell:
    #cdef BASIC_VCELL *thisptr
    def __cinit__(self):
        self.thisptr = new BASIC_VCELL()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr 
