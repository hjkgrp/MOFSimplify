"""
Cython file defining methods for Xyz and Point classes 
declared in geometry.pxd.
"""

__author__ = "Bharat Medasani"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "Dec 09, 2013"


cdef class Xyz:
    """
    Class to store a point
    """
    
    def __cinit__(self, double x=0.0, double y=0.0, double z=0.0):
        self.thisptr = new XYZ(x,y,z)

    def __init__(self, double x=0.0, double y=0.0, double z=0.0):
        pass

    def __dealloc__(self):
        del self.thisptr

    property x:
        def __get__(self): return self.thisptr.x
        def __set__(self, x_in): self.thisptr.x = x_in 
    property y:
        def __get__(self): return self.thisptr.y
        def __set__(self, y_in): self.thisptr.y = y_in 
    property z:
        def __get__(self): return self.thisptr.z
        def __set__(self, z_in): self.thisptr.z = z_in 

    def scale(self, double factor):
        new_xyz = Xyz()
        self.thisptr.scale(factor, new_xyz.thisptr)
        return new_xyz

cdef class Point:
    """
    Class to store a point
    """
    
    def __cinit__(self, double x=0.0, double y=0.0, double z=0.0):
        self.thisptr = new CPoint(x,y,z)

    def __init__(self, double x=0.0, double y=0.0, double z=0.0):
        pass

    def __dealloc__(self):
        del self.thisptr

    def __repr__(self):
        return "("+str(self.x)+','+str(self.y)+','+str(self.y)+')'

    property x:
        def __get__(self): return self.thisptr.vals[0]
        def __set__(self, x_in): self.thisptr.vals[0] = x_in 
    property y:
        def __get__(self): return self.thisptr.vals[1]
        def __set__(self, y_in): self.thisptr.vals[1] = y_in 
    property z:
        def __get__(self): return self.thisptr.vals[2]
        def __set__(self, z_in): self.thisptr.vals[2] = z_in 

    #def scale(self, double scaling_factor):
    #    return self.thisptr.scale(scaling_factor)
