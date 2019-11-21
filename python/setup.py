# -*- coding: utf-8 -*-

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the SBCK (Statistical    ##
## Bias Correction Kit). This library makes it possible to perform bias         ##
## correction with non parametric methods, and give some metrics between Sparse ##
## Histogram is high dimensions.                                                ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## SBCK (Statistical Bias Correction Kit). Cette librairie permet d'appliquer   ##
## une correction de biais avec des méthodes non paramétriques, et propose      ##
## diverses metrique entre Histograme Sparse en haute dimension.                ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################


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
	"SBCK.metrics"
]


setup(
	name = "SBCK" ,
	description = "Statistical Bias Correction Kit" ,
	version = "0.1.1" ,
	author = "Yoann Robin" ,
	author_email = "yoann.robin.k@gmail.com" ,
	license = "CeCILL-C" ,
	platforms = [ "linux" , "macosx" ] ,
	requires = [ "numpy" , "scipy" , "matplotlib" ],
	ext_modules = ext_modules,
	install_requires = ['pybind11>=2.2'],
	cmdclass = {'build_ext': BuildExt},
	zip_safe = False,
	packages = list_packages,
	package_dir = { "SBCK" : "SBCK" }
)


