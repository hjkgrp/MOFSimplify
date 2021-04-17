from libcpp.string cimport  string
from netstorage cimport ATOM_NETWORK

cdef extern from "../../psd.h":
    cdef void c_calcPoreSizeDistr "calcPoreSizeDistr"(ATOM_NETWORK *atmnet, 
            ATOM_NETWORK *orgAtomnet, bint highAccuracy, double r_probe_chan, 
            double r_probe, int numSamples, bint excludePockets, 
            string histFile, string pointsFile, string nodeRadiiFile,
            string spheresDistFile, bint visualize, bint overlapsCheck)
