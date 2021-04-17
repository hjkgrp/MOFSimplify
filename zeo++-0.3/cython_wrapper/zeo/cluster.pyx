"""
This module provides functions to simplify high accuracy voronoi network
"""
from __future__ import print_function

__author__ = 'Bharat Medasani'
__copyright__ = ""
__version__ = "0.1"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "alpha"
__date__ = "Dec 6, 2013"

import sys
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as inc

from zeo.netstorage cimport ATOM_NETWORK, VORONOI_NETWORK, VOR_NODE
from zeo.netstorage cimport AtomNetwork, VoronoiNetwork
from zeo.high_accuracy import high_accuracy_atmnet
from zeo.geometry cimport XYZ
from zeo.geometry import Xyz

def warning(*objs):
    print ("WARNING", *objs, files=sys.stderr)

def simplify_highaccuracy_vornet(atmnet):
    """
    Generates and simplifies high accuracy voronoi network 
    """
    cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
    simplify_ha_vornet(c_atmnetptr)


def reduced_highaccuracy_vornodes(atmnet):
    """
    Generates simplified hgh accuracy voronoi network
    """
    cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
    cdef vector[XYZ] xyz_vect
    high_accuracy_vornodes_reduction(c_atmnetptr, &xyz_vect)
    # Conver to list of Xyz
    xyz_list = []

    cdef vector[XYZ].iterator it = xyz_vect.begin()
    while it != xyz_vect.end():
        new_xyz = Xyz((deref(it)).x, (deref(it)).y, (deref(it)).z) #Infefficient
        xyz_list.append(new_xyz)
        inc(it)

    return xyz_list


def pruned_highaccuracy_voronoi_network(atmnet, delta=0.5):
    """
    Prunes hgh accuracy voronoi network by removing voronoi
    nodes close to the center of the bigger atoms.
    """
    ha_atmnet = atmnet.copy()
    high_accuracy_atmnet(ha_atmnet, "MED")
    vornet,fcs = ha_atmnet.perform_voronoi_decomposition()
    cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
    cdef ATOM_NETWORK* c_ha_atmnetptr = (<AtomNetwork?>ha_atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornetptr = (<VoronoiNetwork?>vornet).thisptr
    prune_high_accuracy_voronoi_network(c_vornetptr, c_atmnetptr, 
            c_ha_atmnetptr, delta)
    return vornet

def get_nearest_largest_diameter_highaccuracy_vornode( atmnet, delta=0.25):
    """
    Get the reduced high accuracy voronoi network where only nodes that 
    has the largest diameter and within the cutoff distance to the nodes
    of the low accuracy voronoi network are retained. A one-one mapping
    of high accuracy voronoi nodes and low accuracy nodes is obtained.

    Input:
        atmnet: AtomNetwork object
        delta: cutoff (default = 0.25 angstroms)
    Output:
        Reduced voronoi network
    """
    #generate_simplified_highaccuracy_voronoi_network(atmnet)
    ha_vornet = pruned_highaccuracy_voronoi_network(atmnet, delta=0.7)
    #print ''
    #print '**********ONE DECOMPOSITION.************'
    #print ''
    vornet, fcs = atmnet.perform_voronoi_decomposition()
    cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef VORONOI_NETWORK* c_ha_vornet_ptr = (<VoronoiNetwork?>ha_vornet).thisptr
    red_vornet = VoronoiNetwork()
    #print ''
    #print '*********WORKED TILL HERE*********'
    #print ''
    nearest_largest_diameter_ha_vornet(c_ha_vornet_ptr, c_vornet_ptr, 
            c_atmnet_ptr, red_vornet.thisptr, delta)
    return red_vornet

def generate_simplified_highaccuracy_voronoi_network(atmnet,delta=0.6):
    """
    Generate a simplified high accuracy voronoi network. 
    Uses Zeo++ high accuracy network and simplifies it such that only voronoi 
    nodes that belong to different atoms of original atom network are 
    retained. There can be different no. of voronoi nodes when compared with
    the voronoi nodes obtained with regular tesselation. 
    Input:
        atmnet: AtomNetwork object
    Output:
        Simplified high accuracy voronoi network
    """
    ha_atmnet = atmnet.copy()
    high_accuracy_atmnet(ha_atmnet, "LOW")
    vornet,ecs,fcs = atmnet.perform_voronoi_decomposition()
    ha_vornet,ecs,fcs = ha_atmnet.perform_voronoi_decomposition()
    node_size = vornet.size()
    ha_node_size = ha_vornet.size()
    if node_size == ha_node_size:
        warning('No high accuracy')
    return ha_vornet        # The processing below is eliminated temporarily
    node_size = vornet.size()
    ha_node_size = ha_vornet.size()
    if node_size == ha_node_size:
        # No need for simplification
        return vornet
    
    #print ''
    #print '**********ONE DECOMPOSITION.************'
    #print ''
    cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef ATOM_NETWORK* c_ha_atmnet_ptr = (<AtomNetwork?>ha_atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef VORONOI_NETWORK* c_ha_vornet_ptr = (<VoronoiNetwork?>ha_vornet).thisptr
    red_vornet = VoronoiNetwork()
    #print ''
    #print '*********WORKED TILL HERE*********'
    #print ''
    simplify_high_accuracy_vornet(c_ha_vornet_ptr, c_ha_atmnet_ptr, 
             red_vornet.thisptr)
    #print 'red_vornet size', red_vornet.size()
    further_red_vornet = VoronoiNetwork()
    #print ''
    #print '*********WORKED TILL HERE*********'
    #print ''
    nearest_largest_diameter_ha_vornet(red_vornet.thisptr, c_vornet_ptr, 
            c_atmnet_ptr, further_red_vornet.thisptr, delta)
    #print ''
    #print '*********WORKED TILL HERE TOO*********'
    #print ''
    #print 'further_red_vornet size', further_red_vornet.size()
    return further_red_vornet
    #print '********SIMPLIFIED_VORNET_COMPLETE*******'

def prune_voronoi_network_close_node(atmnet,delta=0.1):
    """
    Generate a pruned high accuracy voronoi network. 
    Uses Zeo++ high accuracy network and simplifies it such that only voronoi 
    nodes that are farther than "delta" are retained. 
    Input:
        atmnet: AtomNetwork object
    Output:
        Simplified high accuracy voronoi network
    """
    ha_atmnet = atmnet.copy()
    high_accuracy_atmnet(ha_atmnet, "MED")
    vornet,ecs,fcs = atmnet.perform_voronoi_decomposition()
    ha_vornet,ecs,fcs = ha_atmnet.perform_voronoi_decomposition()

    node_size = vornet.size()
    ha_node_size = ha_vornet.size()
    print (node_size, ha_node_size)
    if node_size == ha_node_size:
        warning('No high accuracy')
        #return vornet
    
    cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_ha_vornet_ptr = (<VoronoiNetwork?>ha_vornet).thisptr
    red_vornet = VoronoiNetwork()
    #print ''
    #print '*********WORKED TILL HERE*********'
    #print ''
    geometry_pruning(c_ha_vornet_ptr, c_atmnet_ptr, delta, 
            red_vornet.thisptr)
    print (red_vornet.size())
    pruned_vornet = VoronoiNetwork()
    ha_prune_within_atom(red_vornet.thisptr, c_atmnet_ptr,
            delta, pruned_vornet.thisptr)
    print (pruned_vornet.size())
    return pruned_vornet
