#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup
from setuptools.extension import Extension
try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True


includedirs=["../../../Voro++/src"]
libdirs = ["../../../Voro++/src"]

#print use_cython
#use_cython = False
cmd_class = {}

if use_cython:
    ext = '.pyx'
    cmd_class.update({'build_ext':build_ext})
else:
    ext = '.cpp'

netstorage_srcfiles = [
        'zeo/netstorage'+ext, '../networkstorage.cc', 
        '../mindist.cc', '../geometry.cc', '../networkinfo.cc',
        '../networkio.cc', '../grid.cc', '../symbcalc.cc',
        '../string_additions.cc', 
        '../voronoicell.cc', 
        '../networkanalysis.cc', '../graphstorage.cc', '../area_and_volume.cc',
        '../network.cc', '../OMS.cc', '../v_network.cc', '../symmetry.cc',
        '../networkaccessibility.cc', '../channel.cc', '../net.cc', '../ray.cc',
	'../rmsd.cc','../material.cc', '../psd.cc',
        ]
netinfo_srcfiles = ['zeo/netinfo'+ext, '../networkinfo.cc']
netio_srcfiles = [
        'zeo/netio'+ext, '../networkio.cc', 'zeo/netinfo'+ext, 
        #'../networkinfo.cc', 'zeo/string_add.pxd', '../string_additions.cc', 
        '../networkinfo.cc',  '../string_additions.cc', 
        '../grid.cc', '../mindist.cc', '../symbcalc.cc',  '../symmetry.cc',
        '../networkstorage.cc', '../geometry.cc', '../net.cc', '../rmsd.cc',
        ]
graphstorage_srcfiles = ['zeo/graphstorage'+ext, '../graphstorage.cc']
psd_srcfiles = ['zeo/psd'+ext, '../psd.cc']
voronoicell_srcfiles = [
        'zeo/voronoicell'+ext, '../voronoicell.cc', '../geometry.cc',
	'../networkstorage.cc', '../net.cc', '../mindist.cc', 
	'../networkinfo.cc', '../rmsd.cc', '../symmetry.cc', 
	'../string_additions.cc', '../ray.cc', '../channel.cc', 
	'../network.cc', '../OMS.cc', '../area_and_volume.cc', '../networkaccessibility.cc', 
	'../graphstorage.cc', '../networkanalysis.cc', '../v_network.cc',
        ]
channel_srcfiles = ['zeo/channel'+ext, '../channel.cc']
highaccuracy_srcfiles = [
        'zeo/high_accuracy'+ext, '../sphere_approx.cc', '../networkstorage.cc', 
        '../networkinfo.cc', '../mindist.cc', '../geometry.cc', '../net.cc', 
	'../symmetry.cc', '../string_additions.cc', '../ray.cc',
	'../networkaccessibility.cc', '../network.cc', '../networkio.cc',
	'../grid.cc', '../symbcalc.cc', '../voronoicell.cc', '../graphstorage.cc',
	'../channel.cc', '../v_network.cc', '../networkanalysis.cc',
	'../area_and_volume.cc', '../rmsd.cc', '../material.cc', '../psd.cc',
        ]
areavol_srcfiles = [
        'zeo/area_volume'+ext, '../area_and_volume.cc', '../networkinfo.cc', 
        '../networkstorage.cc', '../mindist.cc', '../geometry.cc', 
        '../networkio.cc', '../grid.cc', '../symbcalc.cc',
        '../string_additions.cc', 'zeo/voronoicell'+ext, '../voronoicell.cc', 
        '../networkanalysis.cc', '../graphstorage.cc', '../symmetry.cc', 
        '../network.cc', '../OMS.cc', '../v_network.cc', '../ray.cc', '../rmsd.cc',
        '../networkaccessibility.cc', '../channel.cc', '../net.cc'
        ]
cluster_srcfiles = [
        'zeo/cluster'+ext, '../cluster.cc', '../networkstorage.cc',
        '../networkinfo.cc', '../mindist.cc', '../geometry.cc', 
        '../network.cc', '../OMS.cc', '../voronoicell.cc', '../graphstorage.cc',
        '../networkanalysis.cc', '../channel.cc', '../v_network.cc',
        '../area_and_volume.cc',  '../networkaccessibility.cc', 
        '../string_additions.cc', '../sphere_approx.cc', '../net.cc',
	'../symmetry.cc', '../ray.cc', '../rmsd.cc', '../material.cc', '../psd.cc',
        ]
cycle_srcfiles = [
        'zeo/cycle'+ext, '../cycle.cc', '../networkstorage.cc',
        '../networkinfo.cc', '../mindist.cc', '../geometry.cc', 
        '../network.cc', '../OMS.cc', '../voronoicell.cc', '../graphstorage.cc',
        '../networkanalysis.cc', '../channel.cc', '../v_network.cc',
        '../area_and_volume.cc',  '../networkaccessibility.cc', 
        '../string_additions.cc', '../sphere_approx.cc'
        ]
geometry_srcfiles = ['zeo/geometry'+ext, '../geometry.cc']
setup(
    name = 'zeo',
    version = '0.1',
    description = "Python interface to Zeo++",
    url = "http://www.maciejharanczyk.info/Zeopp/",
    author = "Bharat Medasani, Maciej Haranzcyk",
    author_email = "bkmedasani@lbl.gov",
    license = "",
    #cmdclass = {'build_ext':build_ext},
    cmdclass = cmd_class,
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Cython",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering"
        ],
    ext_modules = [
                   Extension("zeo.netstorage",
                             sources=netstorage_srcfiles, 
                             include_dirs=includedirs,
                             libraries = ['voro++'], 
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.geometry", 
                             sources=geometry_srcfiles,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.netinfo", 
                             sources=netinfo_srcfiles,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.voronoicell", 
                             sources=voronoicell_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'], 
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.netio",
                             sources=netio_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.graphstorage",
                             sources=graphstorage_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.psd", 
                             sources=psd_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.channel", 
                             sources=channel_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.high_accuracy", 
                             sources=highaccuracy_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.area_volume", 
                             sources=areavol_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.cluster", 
                             sources=cluster_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.cycle", 
                             sources=cycle_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                             ]
)
