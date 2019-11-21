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
import scipy.spatial.distance as ssd
from .__tools_cpp import network_simplex


#############
## Classes ##
#############

class OTHist:##{{{
	def __init__( self , c , p ):
		self.p = p
		self.c = c
##}}}

class OTNetworkSimplex:##{{{
	"""
	SBCK.tools.NetworkSimplex
	=========================
	
	Network simplex method to solve optimal transport problem
	
	References
	==========
	Bazaraa, M. S., Jarvis, J. J., and Sherali, H. D.: Linear Programming and Network Flows, 4th edn., John Wiley & Sons, 2009.
	"""
	
	def __init__( self , power = 2 ):##{{{
		"""
		Initialisation of solver
		
		Parameters
		----------
		power : float
			Power of the plan (default = 2)
		"""
		self.C = None
		self.P = None
		self.power = power
	##}}}
	
	def fit( self , mu0 , mu1 , C = None ):##{{{
		"""
		Fit optimal plan from two measures mu0 and mu1. The measures must have the arguments "c" (center of bins) and "p" (probability of bins)
		
		Parameters
		----------
		mu0 : (SBCK.SparseHist)
			Source histogram
		mu1 : (SBCK.SparseHist)
			Target histogram
		"""
		self.C = ssd.cdist( mu0.c , mu1.c )**self.power if C is None else C
		
		self.P,self.state = network_simplex( mu0.p , mu1.p , self.C )
		
	##}}}
	
	def plan(self):##{{{
		"""
		Return plan estimated
		
		Return
		------
		P : np.array
			Plan
		"""
		return self.P
	##}}}
##}}}

class OTSinkhorn:##{{{
	"""
	SBCK.tools.OTSinkhorn
	=====================
	
	Sinkhorn method to solve optimal transport problem
	
	References
	==========
	Sinkhorn Distances: Lightspeed Computation of Optimal Transportation Distances. arXiv, https://arxiv.org/abs/1306.0895
	
	"""
	
	def __init__( self , power = 2 , eps = 0.1 , tol = 1e-6 ):##{{{
		"""
		Initialisation of solver
		
		Parameters
		----------
		power : float
			Power of the plan (default = 2)
		eps   : float
			Regularization parameter
		tol   : float
			Numerical tolerance
		"""
		self.power = power
		self.eps = eps
		self.tol = tol
		self.u = None
		self.v = None
		self.C = None
		self.K = None
		self.P = None
	##}}}
	
	def fit( self , mu0 , mu1 , C = None ):##{{{
		"""
		Fit optimal plan from two measures mu0 and mu1. The measures must have the arguments "c" (center of bins) and "p" (probability of bins)
		
		Parameters
		----------
		mu0 : (SBCK.SparseHist)
			Source histogram
		mu1 : (SBCK.SparseHist)
			Target histogram
		"""
		self.u = np.ones_like( mu0.p )
		self.v = np.ones_like( mu1.p )
		self.C = ssd.cdist( mu0.c , mu1.c )**self.power if C is None else C
		self.K = np.exp( - self.C / self.eps )
		
		err = 1. + self.tol
		
		while err > self.tol:
			## Iteration
			self.u = mu0.p / ( self.K @ self.v )
			self.v = mu1.p / ( self.K.T @ self.u )
			
			## Margins
			hpX = self.proj0()
			hpY = self.proj1()
			
			## Update error
			err = max( np.linalg.norm( mu0.p - hpX ) , np.linalg.norm( mu1.p - hpY ) )
		
		self.P = np.diag(self.u) @ self.K @ np.diag(self.v)
	##}}}
	
	def proj0( self ):##{{{
		"""
		Projection on first margin
		
		Returns
		-------
		hp0 : np.array
			Vector of probability of the first margin (from plan)
		"""
		return self.u * ( self.K @ self.v   )
	##}}}
	
	def proj1( self ):##{{{
		"""
		Projection on second margin
		
		Returns
		-------
		hp1 : np.array
			Vector of probability of the second margin (from plan)
		"""
		return self.v * ( self.K.T @ self.u )
	##}}}
	
	def plan( self ):##{{{
		"""
		Return plan estimated
		
		Return
		------
		P : np.array
			Plan
		"""
		return np.diag(self.u) @ self.K @ np.diag(self.v)
	##}}}
	
##}}}

class OTSinkhornLogDual:##{{{
	"""
	SBCK.tools.OTSinkhornLogDual
	============================
	
	Sinkhorn method to solve optimal transport problem, on dual problem with min max acceleration, more robust.
	
	"""
	
	def __init__( self , power = 2 , eps = 0.1 , tol = 1e-6 ):##{{{
		"""
		Initialisation of solver
		
		Parameters
		----------
		power : float
			Power of the plan (default = 2)
		eps : float
			Regularization parameter
		tol : float
			Numerical tolerance
		"""
		self.power = power
		self.eps = eps
		self.tol = tol
		self.f = None
		self.g = None
		self.C = None
		
		self.P   = None
		self.hp0 = None
		self.hp1 = None
	##}}}
	
	def fit( self , mu0 , mu1 ):##{{{
		"""
		Fit optimal plan from two measures mu0 and mu1. The measures must have the arguments "c" (center of bins) and "p" (probability of bins)
		
		Parameters
		----------
		mu0 : (SBCK.SparseHist)
			Source histogram
		mu1 : (SBCK.SparseHist)
			Target histogram
		"""
		self.f = np.zeros_like( mu0.p )
		self.g = np.zeros_like( mu1.p )
		self.C = ssd.cdist( mu0.c , mu1.c )**self.power
		
		
		err = 1. + self.tol
		nit = 0
		while err > self.tol and nit < 1000:
			## Iteration of f
			mg     = self.g.reshape( (1,-1) ) - self.C
			maxg   = np.max( mg , axis = 1 ).reshape(-1,1)
			self.f = - maxg.ravel() - self.eps * np.log( np.sum( mu1.p.reshape( (1,-1) ) * np.exp( (mg - maxg) / self.eps ) , axis = 1 ) )
			
			## Iteration of g
			mf     = self.f.reshape( (-1,1) ) - self.C
			maxf   = np.max( mf , axis = 0 ).reshape( (1,-1) )
			self.g = - maxf.ravel() - self.eps * np.log( np.sum( mu0.p.reshape( (-1,1) ) * np.exp( ( mf - maxf ) / self.eps ) , axis = 0 ) )
			
			## Find plan
			self.P = mu0.p.reshape( (-1,1) ) * mu1.p.reshape( (1,-1) ) * np.exp( ( self.f.reshape( (-1,1) ) + self.g.reshape( (1,-1) ) - self.C ) / self.eps )
			
			## Find margins
			self.hp0 = np.sum( self.P , axis = 1 )
			self.hp1 = np.sum( self.P , axis = 0 )
			
			## Update error
			err = max( np.linalg.norm( mu0.p - self.hp0 ) , np.linalg.norm( mu1.p - self.hp1 ) )
			nit += 1
	##}}}
	
	def proj0( self ):##{{{
		"""
		Projection on first margin
		
		Returns
		-------
		hp0 : np.array
			Vector of probability of the first margin (from plan)
		"""
		return self.hp0
	##}}}
	
	def proj1( self ):##{{{
		"""
		Projection on second margin
		
		Returns
		-------
		hp1 : np.array
			Vector of probability of the second margin (from plan)
		"""
		return self.hp1
	##}}}
	
	def plan( self ):##{{{
		"""
		Return plan estimated
		
		Return
		------
		P : np.array
			Plan
		"""
		return self.P
	##}}}
	
##}}}


