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

##################################################################################
##################################################################################
##                                                                              ##
## Original authors : Mathieu Vrac and Soulivanh Thao                           ##
## Contact          : mathieu.vrac@lsce.ipsl.fr                                 ##
## Contact          : soulivanh.thao@lsce.ipsl.fr                               ##
##                                                                              ##
## Notes   : SchaakeShuffleRef is the re-implementation of the rank shuffle     ##
##           funtions of R package "R2D2" developped by Mathieu Vrac and        ##
##           Soulivanh Thao, available at                                       ##
##                                                                              ##
##           This code is governed by the CeCILL-C license with the             ##
##           authorization of Mathieu Vrac                                      ##
##                                                                              ##
##################################################################################
##################################################################################

###############
## Libraries ##
###############

import numpy       as np
import scipy.stats as sc


#############
## Classes ##
#############

class SchaakeShuffle:##{{{
	"""
	SBCK.tools.SchaakeShuffle
	=========================
	Match the rank structure of X with them of Y by reordering X. Work in multivariate case, but rank  of each 
	features are reordered independantly.
	"""
	
	def __init__( self , Y0 = None ):##{{{
		"""
		Initialization
		
		Parameters
		----------
		Y0 : np.array[n_samples,n_features] or None
			The target dataset, use fit function if Y0 is None
		"""
		self._Y0 = None
		if Y0 is not None: self.fit(Y0)
	##}}}
	
	def fit( self , Y0 ):##{{{
		"""
		Definine the reference ranks structure
		
		Parameters
		----------
		Y0 : np.array[n_samples,n_features]
			The target dataset
		"""
		self._Y0 = Y0
		if self._Y0.ndim == 1: self._Y0 = self._Y0.reshape(-1,1)
	##}}}
	
	def _predict( self , Y0 , X0 ):##{{{
		X0 = X0.squeeze()
		Y0 = Y0.squeeze()
		
		rank_X0 = sc.rankdata( X0 , method = "ordinal" )
		rank_Y0 = sc.rankdata( Y0 , method = "ordinal" )
		
		arank_X0 = np.argsort(rank_X0)
		Z0 = X0[arank_X0][rank_Y0-1]
		
		return Z0
	##}}}
	
	def predict( self , X0 ):##{{{
		"""
		Apply the rank structure to X0
		
		Parameters
		----------
		X0 : np.array[n_samples,n_features]
			The dataset to reorder
		
		Returns
		-------
		Z0 : np.array[n_samples,n_features]
			Reordered dataset
		"""
		if X0.ndim == 1: X0 = X0.reshape(-1,1)
		
		## If n_samples of X/Y differs, we complet the sequence by drawing uniformly in X/Y to have the same shape.
		if self._Y0.shape[0] < X0.shape[0]:
				YY = np.zeros_like(X0)
				YY[:self._Y0.shape[0],:] = self._Y0
				YY[self._Y0.shape[0]:,:] = self._Y0[np.random.choice( self._Y0.shape[0] , X0.shape[0] - self._Y0.shape[0] , replace = True ),:]
				XX = X0
		elif X0.shape[0] < self._Y0.shape[0]:
			XX = np.zeros_like(self._Y0)
			XX[:X0.shape[0],:] = X0
			XX[X0.shape[0]:,:] = X0[np.random.choice( X0.shape[0] , self._Y0.shape[0] - X0.shape[0] , replace = True ),:]
			YY = self._Y0
		else:
			XX,YY = X0,self._Y0
		
		n_features = X0.shape[1]
		ZZ = np.zeros_like(XX)
		for i in range(n_features):
			ZZ[:,i] = self._predict( YY[:,i] , XX[:,i] )
		
		Z0 = ZZ[:X0.shape[0],:]
		
		return Z0
	##}}}
##}}}

def schaake_shuffle( Y0 , X0 ):##{{{
	"""
	SBCK.tools.schaake_shuffle
	==========================
	Match the rank structure of X0 with them of Y0 by reordering X0. Work in multivariate case, but ranks of each 
	features are reordered independantly.
	
	Note: This function just call the class SBCK.tools.SchaakeShuffle
	
	Parameters
	----------
	X0 : np.array[n_samples,n_features]
		The dataset to reorder
	Y0 : np.array[n_samples,n_features]
		The target dataset
	
	Returns
	-------
	Z0 : np.array[n_samples,n_features]
		Reordered dataset
	
	"""
	ss = SchaakeShuffle(Y0)
	return ss.predict(X0)
##}}}

class SchaakeShuffleRef(SchaakeShuffle):##{{{
	"""
	SBCK.tools.SchaakeShuffleRef
	============================
	Match the rank structure of X with them of Y by reordering X, but fix one features to keep the structure of X.
	"""
	
	def __init__( self , ref , Y0 = None ):##{{{
		"""
		Initialization
		
		Parameters
		----------
		ref : int
			features kept.
		Y0 : np.array[n_samples,n_features] or None
			The target dataset, use fit function if Y0 is None
		"""
		self._ref = ref
		SchaakeShuffle.__init__( self , Y0 )
	##}}}
	
	def fit( self , Y0 ):##{{{
		"""
		Definine the reference ranks structure
		
		Parameters
		----------
		Y0 : np.array[n_samples,n_features]
			The target dataset
		"""
		SchaakeShuffle.fit( self , Y0 )
	##}}}
	
	def predict( self , X0 ):##{{{
		"""
		Apply the rank structure to X0
		
		Parameters
		----------
		X0 : np.array[n_samples,n_features]
			The dataset to reorder
		
		Returns
		-------
		Z0 : np.array[n_samples,n_features]
			Reordered dataset
		"""
		if X0.ndim == 1: X0 = X0.reshape(-1,1)
		Z0 = SchaakeShuffle.predict( self , X0 )
		
		rank_ref_X0  = sc.rankdata( X0[:,self._ref] , method = "ordinal" )
		rank_ref_Z0  = sc.rankdata( Z0[:,self._ref] , method = "ordinal" )
		arank_ref_Z0 = np.argsort(rank_ref_Z0)
		Z0 = Z0[arank_ref_Z0,:][rank_ref_X0-1,:]
		
		return Z0
	##}}}
##}}}


