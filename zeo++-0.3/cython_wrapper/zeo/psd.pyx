
from libcpp.string cimport string
from netstorage cimport AtomNetwork

def calc_pore_size_distribution(atmnet,  channel_radius, probe_radius, 
        mc_sampling_no, hist_file, high_accuracy=False, exclude_pockets=False, 
        points_file="", node_radii_file="", sphere_dist_file="", 
        vis_flag=False, overlap_check_flag=False):
    """
    Computes the pore size distribution histogram
    Args:
        atmnet:
            zoe.storage.AtomNetwork
        channel_radius:
            Radius of probe used to determine the accessibility of void space.
        probe_radius:
            Radius of probe used in Monte Carlo (MC) sampling of surface.
        mc_sampling_no:
            No. of MC samples per atom
        hist_file:
           File to store the histogram
        high_accuracy (Default=False):
            Optional flag to use high accuracy.
        exclude_pockets (Default=True):
            Optional flag to include pockets.
        points_file (Default=None):
            File to store the points. Used in visualization
        node_radii_file (Default=None):
            File to store the node radi. Used in visualizationi
        sphere_dist_file (Default=None):
            Reserved for future use
        vis_flag (Default=False)
            Visualization Flag
        overlap_check_flag (Default=False)
            VisIT Visualization related Flag
    """
    atmnet_copy = (<AtomNetwork?>atmnet).copy()
    c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    c_atmnetcp_ptr = (<AtomNetwork?>atmnet_copy).thisptr 
    cdef string chist_file = hist_file
    cdef string cpnt_file = points_file
    cdef string cnd_file = node_radii_file
    cdef string csph_file = sphere_dist_file
    c_calcPoreSizeDistr (c_atmnetcp_ptr, c_atmnet_ptr, high_accuracy,
              channel_radius,  probe_radius, mc_sampling_no, exclude_pockets,
              chist_file, cpnt_file, cnd_file, csph_file, vis_flag, 
              overlap_check_flag)

