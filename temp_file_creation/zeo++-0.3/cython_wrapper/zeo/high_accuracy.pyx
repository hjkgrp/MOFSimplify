from zeo.netstorage cimport AtomNetwork, ATOM_NETWORK

_accuracy_kw = {
        "OCC","FCC","ACC","AQC","DDH","TIH","ICH","ICC","RIH","S4","S10","S20",
        "S30","S40","S50","S100","S500","S1000","S10000","DEF","HI","MED","LOW"
        }
def high_accuracy_atmnet(atmnet, accuracy_setting="LOW"):
    """
    Increases the accuracy of voronoi decomposition by replacing big
    atoms (spheres) with a number of small spheres.
    *** Modifies atmnet argument in place ***
    Args:
        atmnet:
            zeo.netstorage.AtomNetwork
            Is modified in place.
        accuracy_setting: 
            String specifying the accuracy settings.
            Possible choices are "OCC","FCC","ACC","AQC","DDH",
            "TIH","ICH","ICC","RIH","S4","S10","S20","S30","S40","S50",
            "S100","S500","S1000","S10000","DEF","HI","MED","LOW".
            Default is "DEF".
    """
    if not accuracy_setting in _accuracy_kw:
        raise ValueError("Accuracy setting not understood")
    cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
    if isinstance(accuracy_setting, unicode):
        accuracy_setting = (<unicode>accuracy_setting).encode('utf8')
    cdef string acc_set = accuracy_setting
    setupHighAccuracyAtomNetwork(c_atmnetptr, acc_set)


