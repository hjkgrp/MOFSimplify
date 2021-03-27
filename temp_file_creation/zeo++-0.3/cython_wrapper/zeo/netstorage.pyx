"""
Cython file defining methods for AtomNetwork and VoronoiNetowrk 
declared in netstorage.pxd file. 
"""

__author__ = "Bharat Medasani"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "Dec 12, 2013"

from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as inc

#from zeo.voronoicell cimport VorCell, BasicVCell, VOR_CELL, BASIC_VCELL
cimport zeo.netinfo
from zeo.voronoicell cimport  VOR_CELL, BASIC_VCELL, VOR_FACE
from zeo.geometry cimport CPoint, Point

#STUFF='Hi'

cdef class Atom:
    """
    Class to store the information about atom (or ion) in a structure.
    """
    def __cinit__(self):
        self.thisptr = new ATOM()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property coords:
        def __get__(self):
            coords = list(self.thisptr.x, self.thisptr.y, self.thisptr.z)
            return coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.x = coords[0]
            self.thisptr.y = coords[1]
            self.thisptr.z = coords[2]

    property radius:
        def __get__(self): return self.thisptr.radius
        def __set__(self, radius): 
            print ("This value is not supposed to be modified")
            self.thisptr.radius = radius


cdef class AtomNetwork:
    """
    Class to store and manipulate the input atom network.
    """
    #Cython wrapper for Zeo++ ATOM_NETWORK class.
    #Contains a pointer to ATOM_NETWORK and a flag denoting whether radius
    #for each atomic species is non-zero. 
    def __cinit__(self):
        self.thisptr = new ATOM_NETWORK()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    def copy(self):
        """
        Create a copy of the AtomNetwork instance
        """
        newatmnet = AtomNetwork()
        self.thisptr.copy(newatmnet.thisptr)
        newatmnet.rad_flag = self.rad_flag
        return newatmnet

    #def relative_to_absolute(self, point):
    #    cdef CPoint* cpoint_ptr = (<Point?>point).thisptr
    #    cdef double x = cpoint_ptr.vals[0]
    #    cdef double y = cpoint_ptr.vals[1]
    #    cdef double z = cpoint_ptr.vals[2]
    #    cdef CPoint abs_point = self.thisptr.abc_to_xyz(x,y,z)
    #    return Point(abs_point.vals[0], abs_point.vals[1], 
    #            abs_point.vals[2])

    #def absolute_to_relative(self, point):
    #    cdef CPoint* cpoint_ptr = (<Point?>point).thisptr
    #    cdef CPoint rel_point = self.thisptr.xyz_to_abc(abs_point)
    #    return Point(rel_point.vals[0], rel_point.vals[1], 
    #            rel_point.vals[2])

    @classmethod
    def read_from_CIF(cls, filename, rad_flag=True, rad_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a CIF file.
        Arguments:
            filename: 
                Input CIF file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, Zeo++ default values are used.
        Returns:
            Instance of AtomNetwork
        """
        #Calls Zeo++ readCIFFile function defined in networkio.cc.
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef char* c_rad_file = rad_file
        if rad_flag:
            if not rad_file:
                zeo.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                c_rad_file = rad_file
                zeo.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readCIFFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_ARC(cls, filename, rad_flag=True, rad_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a ARC file.
        Arguments:
            filename: 
                Input ARC file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readARCFile function defined in networkio.cc.
        cdef char* c_rad_file = rad_file
        if rad_flag:
            if not rad_file:
                zeo.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                c_rad_file = rad_file
                zeo.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readARCFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_CSSR(cls, filename, rad_flag=True, rad_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a CSSR file.
        Arguments:
            filename: 
                Input CSSR file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readCSSRFile function defined in networkio.cc.
        cdef char* c_rad_file
        print (rad_flag, rad_file)
        if rad_flag:
            #if not rad_file:
            zeo.netinfo.zeo_initializeRadTable()
            if rad_file:       # rad_file is defined
                c_rad_file = rad_file
                zeo.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readCSSRFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_V1(cls, filename, rad_flag=True, rad_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a V1 file.
        Arguments:
            filename: 
                Input V1 file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readV1File function defined in networkio.cc.
        cdef char* c_rad_file = rad_file
        if rad_flag:
            if not rad_file:
                zeo.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                zeo.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readV1File(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    def write_to_CSSR(self, filename):
        """
        Writes the atom data in AtomNetwork to a CSSR file.
        Arguments:
            filename: 
                Output CSSR file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToCSSR function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToCSSR(c_filename, self.thisptr):
            raise IOError

    def write_to_CIF(self, filename):
        """
        Writes the atom data in AtomNetwork to a CIF file.
        Arguments:
            filename: 
                Output CIF file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToCIF function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToCIF(c_filename, self.thisptr):
            raise IOError

    def write_to_V1(self, filename):
        """
        Writes the atom data in AtomNetwork to a V1 file.
        Arguments:
            filename: 
                Output V1 file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToV1 function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToV1(c_filename, self.thisptr):
            raise IOError

    def write_to_XYZ(self, filename, supercell_flag, 
                     is_duplicate_perimeter_atoms):
        """
        Writes the atom data in AtomNetwork to an XYZ file.
        Arguments:
            filename: 
                Output XYZ file name.
            supercell_flag:
                Flag denoting whether to write 2x2x2 supercell.
            is_duplicate_perimeter_atoms:
                Flag denoting whether perimeter atoms need to be replicated.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToXYZ function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToXYZ(c_filename, self.thisptr, supercell_flag, 
                is_duplicate_perimeter_atoms):
            raise IOError

    def write_to_VTK(self, filename):
        """
        Writes the boundary of unit cell within the AtomNetwork to a VTK file.
        Arguments:
            filename: 
                Output VTK file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToVTK function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToVTK(c_filename, self.thisptr):
            raise IOError

    def write_to_MOPAC(self, filename, supercell_flag):
        """
        Writes the atom data in AtomNetwork to a .mop file.
        Arguments:
            filename: 
                Output MOPAC file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef char* c_filename = filename
        if not writeToMOPAC(c_filename, self.thisptr, supercell_flag):
             raise IOError

    def calculate_free_sphere_parameters(self, filename):
        """
        Computes the diameters of the largest included sphere, free sphere 
        and included sphere along free sphere path. 
        Arguments:
            filename:
                Name of file where the diameters are stored.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        vornet, edge_centers, face_centers = self.perform_voronoi_decomposition()
        cdef char* c_fname = filename
        vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        calculateFreeSphereParameters(vornet_ptr, c_fname, False)

    def perform_voronoi_decomposition(self, saveVorCells=True):
        """
        Performs weighted voronoi decomposition of atoms in the AtomNetwork 
        to analyze void space and generate voronoi nodes, edges and faces.
        Arguments:
            saveVorCells (optional): 
                Flag to denote whether to save the VorCells.
                Reserved for future use, so ignore this.
        Returns:
            Instance of VoronoiNetwork
        """
        #Calls Zeo++ performVoronoiDecomp function defined in network.cc.
        vornet = VoronoiNetwork()  
        cdef vector[VOR_CELL] vcells
        cdef vector[BASIC_VCELL] bvcells
        #print self.rad_flag
        if not performVoronoiDecomp(self.rad_flag, self.thisptr, 
                vornet.thisptr, &vcells, saveVorCells, &bvcells):
            raise ValueError # Change it to appropriate error
        cdef int N

        # Get the edge centers
        edge_centers = []
        cdef vector[VOR_EDGE] vedges = vornet.thisptr.edges
        cdef vector[VOR_NODE] vnodes = vornet.thisptr.nodes
        for i in range(vedges.size()):
            edge_orig =  vedges[i].origin
            edge_end =  vedges[i].ending
            o_vnode = vnodes[edge_orig]
            e_vnode = vnodes[edge_end]
            edge_center = (o_vnode.x + e_vnode.x, \
                           o_vnode.y + e_vnode.y, \
                           o_vnode.z + e_vnode.z)
            edge_center = tuple(x/2 for x in edge_center)
            if edge_center not in edge_centers:
                edge_centers.append(edge_center)



        # Get the vorcells and obtain the face centers
        face_centers = []
        cdef vector[VOR_FACE] vfaces
        cdef vector[CPoint] vertices
        cdef CPoint* cpoint_ptr 
        #cdef map[int, int] id_maps
        cdef vector[int] node_ids
        face_node_ids = set()
        for i in range(vcells.size()):
            vfaces = vcells[i].faces
            for j in range(vfaces.size()):
                node_ids = vfaces[j].node_ids
                node_id_list = []
                for k in range(node_ids.size()):
                    node_id_list.append(node_ids[k])
                node_id_set = frozenset(node_id_list)
                if not node_id_set in face_node_ids:
                    face_node_ids.add(node_id_set)
                    centroid = Point()
                    cpoint_ptr = (<Point?>centroid).thisptr
                    vertices = vfaces[j].vertices
                    for k in range(vertices.size()):
                        centroid.x = centroid.x + vertices[k].vals[0]
                        centroid.y = centroid.y + vertices[k].vals[1]
                        centroid.z = centroid.z + vertices[k].vals[2]
                    centroid.x = centroid.x/vertices.size()
                    centroid.y = centroid.y/vertices.size()
                    centroid.z = centroid.z/vertices.size()
                    face_centers.append(centroid)

        # Convert the Zeo++ Point objects in (x,y,z) tuple objects
        fcs = []
        for center in face_centers:
            cntr = (center.x,center.y,center.z)
            fcs.append(cntr)

        #bvcelllist = []
        # define copy methods for BASIC_VCELL and VOR_CELL methods
        #for i in range(vcells.size()):
        #    pass
            #vorcell = VorCell()
            #vorcell.thisptr = &(vcells[i])
            #vorcelllist.append(vcells[i])
        #for i in range(bvcells.size()):
        #    pass
            #basicvcell = BasicVCell()
            #basicvcell.thisptr = &(bvcells[i])
            #bvcelllist.append(bvcells[i])
        return vornet, edge_centers, fcs


cdef class VoronoiNode:
    """
    Class to store the voronoi nodes with coordinates and radius
    """
    def __cinit__(self):
        self.thisptr = new VOR_NODE()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property coords:
        def __get__(self):
            coords = list(self.thisptr.x, self.thisptr.y, self.thisptr.z)
            return coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.x = coords[0]
            self.thisptr.y = coords[1]
            self.thisptr.z = coords[2]

    property radius:
        def __get__(self): return self.thisptr.rad_stat_sphere
        def __set__(self, rad): 
            print ("This value is not supposed to be modified")
            self.thisptr.rad_stat_sphere = rad

cdef class VoronoiNetwork:
    """
    Class to store the Voronoi network generated from Voronoi decomposition
    of atom network.
    """
    #Cython wrapper for Zeo++ VORONOI_NETWORK class.
    #Contains a pointer to ATOM_NETWORK and a flag denoting whether radisu
    #for each atomic species is non-zero. 
    def __cinit__(self):
        self.thisptr = new VORONOI_NETWORK()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    def size(self):
        return self.thisptr.nodes.size()

    def prune(self, radius):
        """
        Removes the edges that do not allow a sphere to pass.
        Arguments:
            radius:
                Radius of the sphere
        Returns:
            Instance of VoronoiNetwork with edges pruned.
        """
        cdef VORONOI_NETWORK newcvornet = self.thisptr.prune(radius)
        newvornet = VoronoiNetwork()
        newvornet.thisptr = &newcvornet
        return newvornet

    def analyze_writeto_XYZ(self, name, double probeRad, atmnet, 
            int shift_x=0, int shift_y=0, int shift_z=0):
        """
        Create diagrams of 1) Voronoi network and 2) accessible Voronoi 
        network, and write the diagrams in VTK files and the Voronoi 
        networks in XYZ files. Useful for visualizing the Voronoi network.
        Args:
            name:
                Name to be used for output files.
            probeRad:
                Radius of the probe.
            atmnet:
                zeo.netstorage.AtomNetwork
            shift_x (default=0):
                Shift the accessible Voronoi network along x-axis
            shift_y (default=0):
                Shift the accessible Voronoi network along y-axis
            shift_z (default=0):
                Shift the accessible Voronoi network along z-axis
        """
        if isinstance(name, unicode):
            name = (<unicode>name).encode('utf8')

        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
        cdef char* cname = name
        visVoro(name, probeRad, shift_x, shift_y, shift_z, self.thisptr, 
                c_atmnetptr)

    def write_to_XYZ(self, filename, double cutoff_radius=0):
        """
        Write only voronoi node information to XYZ file.
        Args:
            filename:
                string
                Name of file to which voronoi node info is written.
            cutoff_radius:
                float
                Threshold radius (default=0)
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef char* c_filename = filename
        if not writeVornetToXYZ(c_filename, self.thisptr, 
                cutoff_radius):
            raise ValueError

    @classmethod
    def perform_voronoi_decomposition(cls, atmnet, saveVorCells=False):
        """
        Performs weighted voronoi decomposition of atoms in the AtomNetwork 
        to analyze void space and generate voronoi nodes, edges and faces.
        Arguments:
            saveVorCells (optional): 
                Flag to denote whether to save the VorCells.
                Reserved for future use, so ignore this.
        Returns:
            Instance of VoronoiNetwork
        """
        #Calls Zeo++ performVoronoiDecomp function defined in network.cc.
        vornet = VoronoiNetwork()  
        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
        cdef vector[VOR_CELL] vcells
        cdef vector[BASIC_VCELL] bvcells
        #print self.rad_flag
        if not performVoronoiDecomp(atmnet.rad_flag, c_atmnetptr, 
                vornet.thisptr, &vcells, saveVorCells, &bvcells):
            raise ValueError # Change it to appropriate error
        #cdef int N
        #vorcelllist = []
        #bvcelllist = []
        # define copy methods for BASIC_VCELL and VOR_CELL methods
        #for i in range(vcells.size()):
        #    pass
            #vorcell = VorCell()
            #vorcell.thisptr = &(vcells[i])
            #vorcelllist.append(vcells[i])
        #for i in range(bvcells.size()):
        #    pass
            #basicvcell = BasicVCell()
            #basicvcell.thisptr = &(bvcells[i])
            #bvcelllist.append(bvcells[i])
        return vornet

def substitute_atoms(atmnet, substituteSeed, radialFlag):
    """
    Attempt to substitute every other Si atom with Al atom.
    AtomNetwork may only consist of Si and O atoms, where each Si atom 
    must be bonded to exactly 4 oxygen atoms and each oxygen atom must 
    be bonded to exactly 2 Si atoms. Raises exception if the substitution
    is not successful. 
    Args:
        atmnet:
            zeo.netstorage.AtomNetwork
        substiuteSeed:
            Boolean flag to specify whether the seeded Si atom is 
            substituted or not. Since only 2 configurations are possible 
            if the structure is consistent, changing this parameter enables 
            generation of all configurations. 
        radialFlag:
            Boolean flag to specify whether atomic sizes are to be used.
    Returns:
        If successful, returns AtomNetwork instance with Si replaced with Al
        and the number of substitutions. 
    """
    cdef int substitutionNo[1]
    atmnet_copy = AtomNetwork()
    c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    if not c_substituteAtoms(c_atmnet_ptr, atmnet_copy.thisptr, substituteSeed,
            substitutionNo, radialFlag):
        raise ValueError
    subNo = substitutionNo[0]
    return atmnet_copy, subNo



