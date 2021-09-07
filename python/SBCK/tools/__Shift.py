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
import itertools as itt

###########
## Class ##
###########

class Shift:
	"""
	SBCK.tools.Shift
	================
	
	Description
	-----------
	Shift class used to to transform a dataset X to [X[0:size-lag],X[1:size-lag+1],...]
	"""
	
	def __init__( self , lag , method = "row" , ref = 0 ):
		"""
		Initialisation of shift class
		
		Parameters
		----------
		lag    : integer
			Time lag of the shift
		method : string
			Inverse method, "row" or "col"
		ref    : integer
			Reference columns / rows to inverse
		"""
		self.lag = lag
		self.ref = ref
		self.method = method
	
	@property
	def ref( self ):
		return self._ref
	
	@ref.setter
	def ref( self , ref ):
		self._ref = ref % ( self.lag + 1 )
	
	@property
	def method( self ):
		return self._method
	
	@method.setter
	def method( self , _method ):
		self._method = _method if _method == "row" else "col"
	
	def transform( self , X ):
		"""
		Transform X to the shifted Xs with lag
		
		Parameters
		----------
		X : np.array
			dataset to shift
		
		Returns
		-------
		Xs: np.array
			dataset shifted
		"""
		if X.ndim == 1: X = X.reshape(-1,1) ## genericity to always have a matrix X.
		n_samples,n_features = X.shape
		Xs = np.zeros( ( n_samples - self.lag , ( self.lag + 1 ) * n_features ) )
		
		for i in range(self.lag+1):
			db = i * n_features
			de = i * n_features + n_features
			tb = i
			te = n_samples - ( self.lag + 1 ) + i + 1
			Xs[:,db:de] = X[tb:te,:]
		return Xs
	
	def _inverse_by_row( self , Xs ):
		n_features = int( Xs.shape[1] / ( self.lag + 1 ))
		n_samples  = Xs.shape[0] + self.lag
		
		Xi = np.zeros( (n_samples,n_features) )
		for r in itt.chain(range(self.lag+1),[self.ref]):
			idx  = np.arange( r , n_samples - self.lag , self.lag )
			Xs0  = Xs[idx[:-1],:-n_features].reshape(-1,n_features) ## Without last index, because the last is also the first of next
			Xs1  = Xs[idx[-1],:].reshape(-1,n_features)
			Xs01 = np.vstack( (Xs0,Xs1) )
			n_samples_01 = Xs01.shape[0]
			Xi[r:(r+n_samples_01),:] = Xs01
		
		return Xi
	
	def _inverse_by_col( self , Xs ):
		n_features = int( Xs.shape[1] / (self.lag + 1) )
		n_samples = Xs.shape[0] + self.lag
		Xu   = np.zeros( (n_samples,n_features) )
		
		for i in itt.chain(range(self.lag+1),[self.ref]):
			db = i * n_features
			de = i * n_features + n_features
			tb = i
			te = n_samples - ( self.lag + 1 ) + i + 1
			Xu[tb:te,:] = Xs[:,db:de]
		return Xu
	
	def inverse( self , Xs , method = None ):
		"""
		Inverse transform
		
		Parameters
		----------
		Xs  : np.array
			dataset to unshift
		
		Returns
		-------
		X: np.array
			dataset unshifted
		"""
		if method is not None: self.method = method
		if self.method == "col":
			return self._inverse_by_col(Xs)
		else:
			return self._inverse_by_row( Xs )



