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

import numpy       as np
import scipy.stats as sc
from .__QM import QM
from .__QDM import QDM
from .metrics.__wasserstein import wasserstein
from .tools.__SlopeStoppingCriteria import SlopeStoppingCriteria


###########
## Class ##
###########

class MBCn:
	"""
	SBCK.MBCn
	=========
	
	Description
	-----------
	MBCn Bias correction method, see [1]
	
	References
	----------
	[1] Cannon, Alex J.: Multivariate quantile mapping bias correction: an N-dimensional probability density function transform for climate model simulations of multiple variables, Climate Dynamics, nb. 1, vol. 50, p. 31-49, 10.1007/s00382-017-3580-6
	"""
	def __init__( self , bc = QDM , metric = wasserstein , stopping_criteria = SlopeStoppingCriteria , stopping_criteria_params = { "minit" : 20 , "maxit" : 100 , "tol" : 1e-3 } , **kwargs ): ##{{{
		"""
		Initialisation of MBCn.
		
		Parameters
		----------
		bc  : Bias correction method
			Non stationary bias correction method, default is QDM
		metric : callable
			Callable between two matrices, used as criteria to dermined when stopped iteration. Default is Wasserstein distance
		stopping_criteria: a class to determine when stop iteration
			See note
		stopping_criteria_params : dict
			Params of stopping_criteria
		kwargs : others named arguments
			Passed to bias correction method
		
		Note
		----
		MBCn method used an alternance of random rotation and quantile mapping to perform the multivariate bias correction.
		At each step, the metric is used in calibration period to determine if correction is close or not of correction.
		The class SlopeStoppingCriteria compute the slope of time series of callable. When the slope is lower than tol,
		or maxit is atteigned, iterations are stopped. At least minit is performed.
		Attributes
		----------
		"""
		self.iter_stop = stopping_criteria(**stopping_criteria_params)
		self.metric = metric
		self.bc = bc
		self.bc_params = kwargs
		self._lbc  = []
	##}}}
	
	@property
	def maxit(self):##{{{
		return self.iter_stop.maxit
	##}}}
	
	@property
	def nit(self):##{{{
		return self.iter_stop.nit
	##}}}
	
	def fit( self , Y0 , X0 , X1 ):##{{{
		"""
		Fit the MBCn model
		
		Parameters
		----------
		Y0 : np.array[ shape = (n_samples,n_features) ]
			Reference dataset during period 0
		X0 : np.array[ shape = (n_samples,n_features) ]
			Biased dataset during period 0
		X1	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset during period 1
		"""
		
		if Y0.ndim == 1: Y0 = Y0.reshape(-1,1)
		if X0.ndim == 1: X0 = X0.reshape(-1,1)
		if X1.ndim == 1: X1 = X1.reshape(-1,1)
		
		n_features = Y0.shape[1]
		self.iter_stop.initialize()
		
		## Generate orthogonal matrices: SO(n_features)
		self.ortho_mat  = sc.special_ortho_group.rvs( n_features , self.maxit )
		
		## Tips for performance, inverse + ortho of next in one pass
		self.tips = np.zeros(self.ortho_mat.shape)
		for i in range(self.maxit-1):
			self.tips[i,:,:] = self.ortho_mat[i+1,:,:] @ np.linalg.inv(self.ortho_mat[i,:,:])
		self.tips[-1,:,:] = np.linalg.inv(self.ortho_mat[-1,:,:])
		
		## Loop
		Z0_o = np.transpose( self.ortho_mat[0,:,:] @ X0.T )
		Z1_o = np.transpose( self.ortho_mat[0,:,:] @ X1.T )
		
		for i in self.iter_stop:
			Y0_o = np.transpose( self.ortho_mat[i,:,:] @ Y0.T )
			
			bc = self.bc(**self.bc_params)
			bc.fit( Y0_o , Z0_o , Z1_o )
			Z1_o,Z0_o = bc.predict(Z1_o,Z0_o)
			
			self.iter_stop.append(self.metric(Z0_o,Y0_o))
			
			self._lbc.append(bc)
			
			Z0_o = np.transpose( self.tips[i,:,:] @ Z0_o.T )
			Z1_o = np.transpose( self.tips[i,:,:] @ Z1_o.T )
		
		Z0 = np.transpose( np.linalg.inv(self.ortho_mat[self.nit,:,:]) @ Z0_o.T )
		Z1 = np.transpose( np.linalg.inv(self.ortho_mat[self.nit,:,:]) @ Z1_o.T )
		
		self.ortho_mat = self.ortho_mat[:self.nit,:,:]
		self.tips = self.tips[:self.nit,:,:]
		self.tips[-1,:,:] = np.linalg.inv(self.ortho_mat[-1,:,:])
		
		bc = self.bc(**self.bc_params)
		bc.fit( Y0 , Z0 , Z1 )
		self._lbc.append(bc)
	##}}}
	
	def predict( self , X1 , X0 = None ):##{{{
		"""
		Perform the bias correction
		Return Z1 if X0 is None, else return a tuple Z1,Z0
		
		Parameters
		----------
		X1  : np.array[ shape = (n_samples,n_features) ]
			Array of value to be corrected in projection period
		X0  : np.array[ shape = (n_samples,n_features) ] or None
			Array of value to be corrected in calibration period
		
		Returns
		-------
		Z1 : np.array[ shape = (n_sample,n_features) ]
			Return an array of correction in projection period
		Z0 : np.array[ shape = (n_sample,n_features) ] or None
			Return an array of correction in calibration period
		"""
		if X0 is None:
			return self._predict_X1(X1)
		else:
			return self._predict_X1_X0(X1,X0)
		
	##}}}
	
	def _predict_X1( self , X1 ):##{{{
		if X1.ndim == 1: X1 = X1.reshape(-1,1)
		
		Z1_o = np.transpose( self.ortho_mat[0,:,:] @ X1.T ) if X1 is not None else None
		
		for i in range(self.nit):
			Z1_o = self._lbc[i].predict(Z1_o)
			Z1_o = np.transpose( self.tips[i,:,:] @ Z1_o.T )
		
		Z1 = self._lbc[-1].predict(Z1_o)
		return Z1
	##}}}
	
	def _predict_X1_X0( self , X1 , X0 ):##{{{
		if X0.ndim == 1: X0 = X0.reshape(-1,1)
		if X1.ndim == 1: X1 = X1.reshape(-1,1)
		
		Z0_o = np.transpose( self.ortho_mat[0,:,:] @ X0.T )
		Z1_o = np.transpose( self.ortho_mat[0,:,:] @ X1.T )
		
		for i in range(self.nit):
			Z1_o,Z0_o = self._lbc[i].predict(Z1_o,Z0_o)
			Z0_o = np.transpose( self.tips[i,:,:] @ Z0_o.T )
			Z1_o = np.transpose( self.tips[i,:,:] @ Z1_o.T )
		
		return self._lbc[-1].predict(Z1_o,Z0_o)
	##}}}


