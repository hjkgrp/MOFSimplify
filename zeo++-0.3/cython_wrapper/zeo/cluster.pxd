from libcpp.vector cimport vector 

from zeo.netstorage cimport ATOM_NETWORK, VORONOI_NETWORK
from zeo.geometry cimport XYZ

"""
cutoff default value is 0.25
"""
cdef extern from "../../cluster.h":
    cdef void simplify_ha_vornet(ATOM_NETWORK *atmnt)

    cdef void high_accuracy_vornodes_reduction(ATOM_NETWORK*, vector[XYZ]*)

    cdef void prune_high_accuracy_voronoi_network(VORONOI_NETWORK* vornet, 
            ATOM_NETWORK* atmnet, ATOM_NETWORK* ha_atmnet, double delta)

    cdef void nearest_largest_diameter_ha_vornet(VORONOI_NETWORK* ha_vornet,
            VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet, 
            VORONOI_NETWORK* red_vornet, float cutoff)

    cdef void simplify_high_accuracy_vornet(VORONOI_NETWORK* ha_vornet,
            ATOM_NETWORK* ha_atmnet, VORONOI_NETWORK* red_vornet)

    cdef void geometry_pruning(VORONOI_NETWORK* ha_vornet, 
            ATOM_NETWORK* ha_atmnet, float cutoff, VORONOI_NETWORK* red_vornet)

    cdef void ha_prune_within_atom(VORONOI_NETWORK* ha_vornet, 
            ATOM_NETWORK* atmnet, float cutoff, VORONOI_NETWORK* red_vornet)
