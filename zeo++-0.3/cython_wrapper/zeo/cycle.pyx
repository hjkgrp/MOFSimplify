"""
Cython file defining methods for computing the centroids of cycles
"""

from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as inc
from netstorage cimport VoronoiNetwork, VORONOI_NETWORK
from netstorage cimport AtomNetwork, ATOM_NETWORK
from geometry cimport  Xyz, XYZ


def compute_centroid_4cycles(vornet):
    """
    Computes the centroid of the 4 corners of quadrilateral voronoi face
    Args:
        vornet:
            zeo.storage.VoronoiNetwork
    Returns:
        List of centroids in [(x1,y1,z1),(x2,y2,z2),...] format
    """

    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef vector[CYCLE] cycles
    cdef vector[int] ids

    if not compute_4cycle(c_vornet_ptr, &cycles, False, 1):
        raise ValueError

    centroid_list = []
    cdef vector[CYCLE].iterator it = cycles.begin()
    cdef vector[int].iterator iit
    while it != cycles.end():
        new_xyz = Xyz()
        centroid(&(deref(it)), new_xyz.thisptr, &ids)
        iit = ids.begin()
        #print ids.size()
        id_set = set()
        while iit != ids.end():
            id_set.add(deref(iit))
            inc(iit)
        
        centroid_list.append({'ids':id_set, 'coords':new_xyz})
        inc(it)

    return centroid_list

def compute_face_centers(atmnet):
    """
    Compute the face centers of the voronoi network 
    """
    cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef vector[XYZ] points
    face_center(c_atmnet_ptr, &points)


