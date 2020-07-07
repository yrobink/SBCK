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

import numpy as np
import scipy.stats as sc
import scipy.interpolate as sci


#############
## Classes ##
#############

class MonotoneInverse:##{{{
	
	def __init__( self , xminmax , yminmax , transform ):##{{{
		self.xmin  = xminmax[0]
		self.xmax  = xminmax[1]
		self.ymin  = yminmax[0]
		self.ymax  = yminmax[1]
		delta = 0.05 * (self.xmax - self.xmin)
		nstepmin,nstepmax = 0,0
		while transform(self.xmin) > self.ymin:
			self.xmin -= delta
			nstepmin += 1
		while transform(self.xmax) < self.ymax:
			self.xmax += delta
			nstepmax += 1
		self.nstep = 100 + max(nstepmin,nstepmax)
		x = np.linspace(self.xmin,self.xmax,self.nstep)
		y = transform(x)
		self._inverse = sci.interp1d( y , x )
	##}}}
	
	def __call__( self , y ):##{{{
		return self._inverse(y)
	##}}}

##}}}

class rv_histogram(sc.rv_histogram):##{{{
	"""
	SBCK.tools.rv_histogram
	=======================
	Wrapper on scipy.stats.rv_histogram adding a fit method.
	"""
	def __init__( self , *args , **kwargs ):##{{{
		sc.rv_histogram.__init__( self , *args , **kwargs )
	##}}}
	
	def fit( X , bins = 100 ):##{{{
		return (np.histogram( X , bins = bins ),)
	##}}}
	
##}}}

class rv_ratio_histogram(sc.rv_histogram):##{{{
	"""
	SBCK.tools.rv_ratio_histogram
	=============================
	Extension of SBCK.tools.rv_histogram taking into account of a "ratio" part, i.e., instead of fitting:
	P( X < x )
	We fit separatly the frequency of 0 and:
	P( X < x | X > 0 )
	"""
	def __init__( self , *args , **kwargs ):##{{{
		eargs = ()
		if len(args) > 0:
			eargs = (args[0],)
		sc.rv_histogram.__init__( self , *eargs , **kwargs )
		self.p0 = 0
		if len(args) > 1:
			self.p0 = args[1]
	##}}}
	
	def fit( X , bins = 100 ):##{{{
		Xp = X[X>0]
		p0 = np.sum(np.logical_not(X>0)) / X.size
		return (np.histogram( Xp , bins = bins ),p0)
	##}}}
	
	def cdf( self , x ):##{{{
		cdf = np.zeros_like(x)
		idxp = x > 0
		idx0 = np.logical_not(x>0)
		cdf[idxp] = (1-self.p0) * sc.rv_histogram.cdf( self , x[idxp] ) + self.p0
		cdf[idx0] = self.p0 / 2
		return cdf
	##}}}
	
	def ppf( self , p ):##{{{
		idxp = p > self.p0
		idx0 = np.logical_not(p > self.p0 )
		ppf = np.zeros_like(p)
		ppf[idxp] = sc.rv_histogram.ppf( self , (p[idxp] - self.p0) / (1-self.p0) )
		ppf[idx0] = 0
		return ppf
	##}}}
	
	def sf( self , x ):##{{{
		return 1 - self.cdf(x)
	##}}}
	
	def isf( self , p ):##{{{
		return self.ppf( 1 - p )
	##}}}

##}}}

class rv_density:##{{{
	
	def __init__( self , *args , **kwargs ):##{{{
		self._kernel = None
		if kwargs.get("X") is not None:
			X = kwargs.get("X")
			self._kernel = sc.gaussian_kde( X.squeeze() , bw_method = kwargs.get("bw_method") )
			self._init_icdf( [X.min(),X.max()] )
		elif len(args) > 0:
			self._kernel = args[0]
			self._init_icdf( [args[1],args[2]] )
	##}}}
	
	def rvs( self , size ):##{{{
		p = np.random.uniform( size = size )
		return self.icdf(p)
	##}}}
	
	def fit( X , bw_method = None ):##{{{
		kernel = sc.gaussian_kde( X , bw_method = bw_method )
		return (kernel,X.min(),X.max())
	##}}}
	
	def pdf( self , x ):##{{{
		return self._kernel.pdf(x)
	##}}}
	
	def cdf( self , x ):##{{{
		x = np.array([x]).squeeze().reshape(-1,1)
		cdf = np.apply_along_axis( lambda z: self._kernel.integrate_box_1d( -np.Inf , z ) , 1 , x )
		cdf[cdf < 0] = 0
		cdf[cdf > 1] = 1
		return cdf.squeeze()
	##}}}
	
	def sf( self , x ):##{{{
		return 1 - self.cdf(x)
	##}}}
	
	def ppf( self , q ):##{{{
		return self.icdf(q)
	##}}}
	
	def icdf( self , q ):##{{{
		return self._icdf_fct(q)
	##}}}
	
	def isf( self , q ):##{{{
		return self.icdf(1-q)
	##}}}
	
	def _init_icdf( self , xminmax ):##{{{
		self._icdf_fct = MonotoneInverse( xminmax , [0,1] , self.cdf )
	##}}}

##}}}

class rv_mixture:##{{{
	
	def __init__( self , l_dist , weights = None ):##{{{
		self._l_dist  = l_dist
		self._n_dist  = len(l_dist)
		self._weights = np.array([weights]).squeeze() if weights is not None else np.ones(self._n_dist)
		self._weights /= self._weights.sum()
		self._init_icdf()
	##}}}
	
	def rvs( self , size ):##{{{
		out = np.zeros(size)
		ib,ie = 0,int(self._weights[0]*size)
		for i in range(self._n_dist-1):
			out[ib:ie] = self._l_dist[i].rvs( size = ie - ib )
			next_size = int(self._weights[i+1]*size)
			ib,ie = ie,min(ie+next_size,size)
		out[ib:] = self._l_dist[-1].rvs( size = size - ib )
		
		return out[np.random.choice(size,size,replace = False)]
	##}}}
	
	def pdf( self , x ):##{{{
		x = np.array([x]).reshape(-1,1)
		dens = np.zeros_like(x)
		for i in range(self._n_dist):
			dens += self._l_dist[i].pdf(x) * self._weights[i]
		return dens
	##}}}
	
	def cdf( self , x ):##{{{
		x = np.array([x]).reshape(-1,1)
		cdf = np.zeros_like(x)
		for i in range(self._n_dist):
			cdf += self._l_dist[i].cdf(x) * self._weights[i]
		return cdf.squeeze()
	##}}}
	
	def sf( self , x ):##{{{
		return 1 - self.cdf(x)
	##}}}
	
	def ppf( self , q ):##{{{
		return self.icdf(q)
	##}}}
	
	def icdf( self , q ):##{{{
		q = np.array([q]).reshape(-1,1)
		return self._icdf_fct(q)
	##}}}
	
	def _init_icdf(self):##{{{
		rvs = self.rvs(10000)
		self._icdf_fct = MonotoneInverse( [rvs.min(),rvs.max()] , [0,1] , self.cdf )
	##}}}
	
	def isf( self , q ):##{{{
		return self.icdf(1-q)
	##}}}
	
##}}}


