# distutils: language = c++
# distutils: sources = ../graphstorage.cc

from netstorage cimport VoronoiNetwork
cdef class DijkstraNetwork:
    """
    Python wrapper class to Zeo++ Djikstra Network
    """
    #cdef DIJKSTRA_NETWORK* thisptr
    def __cinit__(self):
        self.thisptr = new DIJKSTRA_NETWORK()
    @classmethod
    def from_VoronoiNetwork(vornet):
        """
        Build Dijkstra Net from input Voronoi Net
        """
        dijkstranet = DijkstraNetwork()
        c_vornet = (<VoronoiNetwork?>vornet).thisptr
        buildDijkstraNetwork(c_vornet, dijkstranet.thisptr)
        return dijkstranet 
    def __dealloc__(self):
        del self.thisptr

