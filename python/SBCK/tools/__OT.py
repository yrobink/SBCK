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


