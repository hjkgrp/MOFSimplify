from libcpp.string cimport string

from zeo.netstorage cimport AtomNetwork, ATOM_NETWORK
import zeo.high_accuracy  

def volume(atmnet, channel_radius, probe_radius, 
        mc_sampling_no, high_accuracy=False, high_accuracy_atmnet=None, 
        exclude_pockets=True, low_dist_range=-1, high_dist_range=-1):
    """
    Calculates the volume of channels and pockets in a given strucutre.
    Args:
        atmnet:
            zoe.storage.AtomNetwork
        channel_radius:
            Radius of probe used to determine the accessibility of void space.
        probe_radius:
            Radius of probe used in Monte Carlo (MC) sampling of surface.
        mc_sampling_no:
            No. of MC samples per atom
        high_accuracy (Default=False):
            Optional flag to use high accuracy.
        high_accuracy_atmnet (Default=None):
            zeo.netstorage.AtomNetwork
            Optional high accuracy AtomNetwork. If not given and high_accuracy
            flag is set to True, then it is computed and returned.
        exclude_pockets (Default=True):
            Optional flag to include pockets.
        low_dist_range(Default=-1):
            Use if you know the C++ Zeo++ code.
        high_dist_range(Default=-1):
            Use if you know the C++ Zeo++ code.
    Returns:
        1) string containing channel and pocket volumes
        2) if high_accuracy=True and no input high_accuracy_atmnet is given,
           returns high_accuracy_atmnet for future use.
    """
    if high_accuracy and not high_accuracy_atmnet:
        high_accuracy_atmnet = atmnet.copy()
        zeo.high_accuracy.high_accuracy_atmnet(high_accuracy_atmnet)
        ret_high_acc_atmnet = True
    else:
        ret_high_acc_atmnet = False

    if high_accuracy_atmnet and not high_accuracy:
        high_accuracy = True

    cdef ATOM_NETWORK* c_org_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef ATOM_NETWORK* c_atmnet_ptr
    if high_accuracy_atmnet:
        c_atmnet_ptr = (<AtomNetwork?>high_accuracy_atmnet).thisptr
    else:
        tmp_atmnet = atmnet.copy()
        c_atmnet_ptr = (<AtomNetwork?>tmp_atmnet).thisptr

    vol_str = calcAV(c_atmnet_ptr, c_org_atmnet_ptr, high_accuracy, 
            channel_radius, probe_radius, mc_sampling_no, exclude_pockets,
            low_dist_range, high_dist_range)
    #print vol_str
    if ret_high_acc_atmnet:
        return vol_str, high_accuracy_atmnet
    else:
        return vol_str

    #lines = vol_str.split('\n')
    #for line in lines:
    #    if "Number_of_pockets" in line
    #        fields = line.split(" ")
    #        print fields[1], fields[3]

def surface_area(atmnet, channel_radius, probe_radius, 
        mc_sampling_no, high_accuracy=False, high_accuracy_atmnet=None, 
        exclude_pockets=True, extended_output=False):

    """
    Calculates the surface area of channels and pockets in a given strucutre.
    Args:
        atmnet:
            zoe.storage.AtomNetwork
        channel_radius:
            Radius of probe used to determine the accessibility of void space.
        probe_radius:
            Radius of probe used in Monte Carlo (MC) sampling of surface.
        mc_sampling_no:
            No. of MC samples per atom
        high_accuracy (Default=False):
            Optional flag to use high accuracy.
        high_accuracy_atmnet (Default=None):
            zeo.netstorage.AtomNetwork
            Optional high accuracy AtomNetwork. If not given and high_accuracy
            flag is set to True, then it is computed and returned.
        exclude_pockets (Default=True):
            Optional flag to include pockets.
        low_dist_range(Default=-1):
            Use if you know the C++ Zeo++ code.
        high_dist_range(Default=-1):
            Use if you know the C++ Zeo++ code.
    Returns:
        1) string containing channel and pocket surface area
        2) if high_accuracy=True and no input high_accuracy_atmnet is given,
           returns high_accuracy_atmnet for future use.
    """
    if high_accuracy and not high_accuracy_atmnet:
        high_accuracy_atmnet = atmnet.copy()
        zeo.high_accuracy.high_accuracy_atmnet(high_accuracy_atmnet)
        ret_high_acc_atmnet = True
    else:
        ret_high_acc_atmnet = False

    if high_accuracy_atmnet and not high_accuracy:
        high_accuracy = True

    cdef ATOM_NETWORK* c_org_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef ATOM_NETWORK* c_atmnet_ptr
    if high_accuracy_atmnet:
        c_atmnet_ptr = (<AtomNetwork?>high_accuracy_atmnet).thisptr
    else:
        tmp_atmnet = atmnet.copy()
        c_atmnet_ptr = (<AtomNetwork?>tmp_atmnet).thisptr

    sa_str = calcASA(c_atmnet_ptr, c_org_atmnet_ptr, high_accuracy,
            channel_radius, probe_radius, mc_sampling_no, exclude_pockets,
            extended_output)
    if ret_high_acc_atmnet:
        return sa_str, high_accuracy_atmnet
    else:
        return sa_str
