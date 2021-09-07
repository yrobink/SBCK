# -*- coding: utf-8 -*-

## Copyright(c) 2021 Yoann Robin
## 
## This file is part of SBCK.
## 
## SBCK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## SBCK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with SBCK.  If not, see <https://www.gnu.org/licenses/>.


###############
## Libraries ##
###############

import sys,os
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools


#####################
## User Eigen path ##
#####################

eigen_usr_include = ""

i_eigen = -1
for i,arg in enumerate(sys.argv):
	if arg[:5] == "eigen":
		eigen_usr_include = arg[6:]
		i_eigen = i

if i_eigen > -1:
	del sys.argv[i_eigen]


################################################################
## Some class and function to compile with Eigen and pybind11 ##
################################################################

class get_pybind_include(object):##{{{
	"""Helper class to determine the pybind11 include path
	The purpose of this class is to postpone importing pybind11
	until it is actually installed, so that the ``get_include()``
	method can be invoked. """
	
	def __init__(self, user=False):
		self.user = user
	
	def __str__(self):
		import pybind11
		return pybind11.get_include(self.user)
##}}}

def get_eigen_include( propose_path = "" ):##{{{
	
	possible_path = [ propose_path , "/usr/include/" , "/usr/local/include/" ]
	if os.environ.get("HOME") is not None:
		possible_path.append( os.path.join( os.environ["HOME"] , ".local/include" ) )
	
	for path in possible_path:
		
		
		eigen_include = os.path.join( path , "Eigen" )
		if os.path.isdir( eigen_include ):
			return path
		
		eigen_include = os.path.join( path , "eigen3" , "Eigen" )
		if os.path.isdir( eigen_include ):
			return os.path.join( path , "eigen3" )
	
	return ""
##}}}

def has_flag(compiler, flagname):##{{{
	"""Return a boolean indicating whether a flag name is supported on
	the specified compiler.
	"""
	import tempfile
	with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
		f.write('int main (int argc, char **argv) { return 0; }')
		try:
			compiler.compile([f.name], extra_postargs=[flagname])
		except setuptools.distutils.errors.CompileError:
			return False
	return True
##}}}

def cpp_flag(compiler):##{{{
	"""Return the -std=c++[11/14] compiler flag.
	The c++14 is prefered over c++11 (when it is available).
	"""
	if has_flag(compiler, '-std=c++14'):
		return '-std=c++14'
	elif has_flag(compiler, '-std=c++11'):
		return '-std=c++11'
	else:
		raise RuntimeError( 'Unsupported compiler -- at least C++11 support is needed!' )
##}}}

class BuildExt(build_ext):##{{{
	"""A custom build extension for adding compiler-specific options."""
	c_opts = {
		'msvc': ['/EHsc'],
		'unix': [],
	}
	
	if sys.platform == 'darwin':
		c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
	
	def build_extensions(self):
		ct = self.compiler.compiler_type
		opts = self.c_opts.get(ct, [])
		opts.append( "-O3" )
		if ct == 'unix':
			opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
			opts.append(cpp_flag(self.compiler))
			if has_flag(self.compiler, '-fvisibility=hidden'):
				opts.append('-fvisibility=hidden')
		elif ct == 'msvc':
			opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
		for ext in self.extensions:
			ext.extra_compile_args = opts
		build_ext.build_extensions(self)
##}}}


##########################
## Extension to compile ##
##########################

ext_modules = [
	Extension(
		'SBCK.tools.__tools_cpp',
		['SBCK/tools/src/tools.cpp'],
		include_dirs=[
			# Path to pybind11 headers
			get_eigen_include(eigen_usr_include),
			get_pybind_include(),
			get_pybind_include(user=True)
		],
		language='c++',
		depends = [
			"SBCK/tools/src/SparseHist.hpp"
			"SBCK/tools/src/NetworkSimplex.hpp"
			"SBCK/tools/src/NetworkSimplexLemon.hpp"
			]
	),
]


#################
## Compilation ##
#################

list_packages = [
	"SBCK",
	"SBCK.tools",
	"SBCK.metrics",
	"SBCK.datasets"
]


setup(
	name = "SBCK" ,
	description = "Statistical Bias Correction Kit" ,
	version = "0.4.0" ,
	author = "Yoann Robin" ,
	author_email = "yoann.robin.k@gmail.com" ,
	license = "GNU-GPL3" ,
	platforms = [ "linux" , "macosx" ] ,
	requires = [ "numpy" , "scipy" , "matplotlib" ],
	ext_modules = ext_modules,
	install_requires = ['pybind11>=2.2'],
	cmdclass = {'build_ext': BuildExt},
	zip_safe = False,
	packages = list_packages,
	package_dir = { "SBCK" : "SBCK" }
)


