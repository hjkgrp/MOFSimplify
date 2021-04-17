#!/usr/bin/env python

import os
from setuptools import setup
from setuptools.extension import Extension
#from distutils.extension import Extension
#from distutils.core import setup
from Cython.Distutils import build_ext
from Cython.Build import cythonize

includedirs=["../../../../Voro++/voro/trunk/src", ".."]
libdirs = ["../../../../Voro++/voro/trunk/src", ".."]
runtimedir = os.path.realpath("..")
netstorage_srcfiles = ['zeo/netstorage.pyx' ]
netinfo_srcfiles = ['zeo/netinfo.pyx']
netio_srcfiles = ['zeo/netio.pyx',  'zeo/netinfo.pyx']
graphstorage_srcfiles = ['zeo/graphstorage.pyx']
psd_srcfiles = ['zeo/psd.pyx']
voronoicell_srcfiles = ['zeo/voronoicell.pyx']
channel_srcfiles = ['zeo/channel.pyx']
highaccuracy_srcfiles = ['zeo/high_accuracy.pyx']
areavol_srcfiles = ['zeo/area_volume.pyx']
cluster_srcfiles = ['zeo/cluster.pyx' ]
geometry_srcfiles = ['zeo/geometry.pyx']
cycle_srcfiles = ['zeo/cycle.pyx']

setup(
    name = 'zeo',
    version = '0.1',
    description = "Python interface to Zeo++",
    url = "http://www.maciejharanczyk.info/Zeopp/",
    author = "Bharat Medasani, Maciej Haranzcyk",
    author_email = "bkmedasani@lbl.gov",
    license = "",
    cmdclass = {'build_ext':build_ext},
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Cython",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering"
        ],
    ext_modules = [Extension("zeo.voronoicell", 
                             sources=voronoicell_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'], 
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             #extra_compiler_args = ['-fPIC -Wl,-rpath-link=/home/mbkumar/lbnl/Zeo++/zeo/trunk/'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.netstorage",
                             sources=netstorage_srcfiles, 
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'], 
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.netinfo", 
                             sources=netinfo_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'], 
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.netio",
                             sources=netio_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'],
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.graphstorage",
                             sources=graphstorage_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'], 
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.psd", 
                             sources=psd_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'],
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.channel", 
                             sources=channel_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'], 
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.high_accuracy", 
                             sources=highaccuracy_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'],
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.area_volume", 
                             sources=areavol_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'],
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.cluster", 
                             sources=cluster_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'],
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.geometry", 
                             sources=geometry_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'],
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                   Extension("zeo.cycle", 
                             sources=cycle_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++', 'zeo++'],
                             library_dirs = libdirs,
                             #extra_compiler_args = ['-fPIC'],
                             extra_compiler_args = ['-fPIC -Wl,-rpath='+runtimedir],
                             runtime_library_dirs=[runtimedir],
                             language = 'c++'),
                             ]
)
