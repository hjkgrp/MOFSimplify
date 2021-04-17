# distutils: language = c++
# distutils: sources = ../networkinfo.cc
"""
Wrapper functions to Zeo++ atomic definitons and related functions
"""

#Python definitions for the cdefinitions in .pxd file
def initializeRadTable():
    """
    Populate the atomic radius table with Zeo++ default values
    """
    zeo_initializeRadTable()

def initializeCovRadTable():
    """
    Populate the covalent tradius table with Zeo++ default values
    """
    zeo_initializeCovRadTable()

def initializeMassTable():
    """
    Populate the atomic mass table with Zeo++ default values
    """
    zeo_initializeMassTable()

def initializeAtomCharacterTable():
    """
    Populate the Atom symbol table with Zeo++ default values
    """
    zeo_initializeAtomCharacterTable()

def initializeAtomicNumberTable():
    """
    Populate the atomic number table with Zeo++ default values
    """
    zeo_initializeAtomicNumberTable()

def readRadTable(filename):
    """
    Read atomic radii values from input file and replace the default values
    """
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    zeo_readRadTable(c_filename)

def readMassTable(filename):
    """
    Read atomic mass values from input file and replace the default values
    """
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    zeo_readMassTable(c_filename)

def lookupRadius(element):
    """"
    Args:
        element:
            Element name in conventional shorthand 
            Ex: Al for aluminum 
                Si for silicon 
    Returns:
        radius of the input element
    """
    radius = zeo_lookupRadius(element, True)
    return radius
    
def lookupCovRadius(element):
    return zeo_lookupCovRadius(element)

def lookupMass(element):
    return zeo_lookupMass(element)

def lookupAtomicNumber(element):
    return zeo_lookupAtomicNumber(element)

def isMetal(element):
    return zeo_isMetal(element)
