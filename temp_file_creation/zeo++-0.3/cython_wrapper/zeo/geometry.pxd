"""
Cython declarations file for Zeo++ classes in geometry module.
Declares Zeo++ XYZ and Point Classes and the associated Python wrappers.
"""
__author__ = "Bharat Medasani"
__date__ = "2013-12-09"

cdef extern from "../../geometry.h":
    cdef cppclass XYZ:
        XYZ(double, double, double) except +
        double x, y, z
        void scale (const double sc, XYZ*) 

    cdef cppclass CPoint "Point":
        CPoint(double, double, double) except +
        double vals[3]
        CPoint scale (const double) const
        CPoint operator-(CPoint) 
        CPoint operator+(CPoint)
        CPoint operator*(CPoint)



cdef class Xyz:
    """
    Cython wrapper declaration for Zeo++ XYZ class defined in geometry.h
    Contains a pointer to XYZ.
    """
    cdef XYZ* thisptr


cdef class Point:
    """
    Cython wrapper declaration for Zeo++ Point class defined in geometry.h
    Contains a pointer to c_point.
    """
    cdef CPoint* thisptr

