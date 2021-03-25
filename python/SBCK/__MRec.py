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
from .tools.__rv_extend import rv_histogram
from .__QM import QM


###########
## Class ##
###########


class MRec:
	"""
	SBCK.MRec
	=========
	
	Description
	-----------
	MRec Bias correction method, see [1]
	
	References
	----------
	[1] Bárdossy, A. and Pegram, G.: Multiscale spatial recorrelation of RCM precipitation to produce unbiased climate change scenarios over large areas and small, Water Resources Research, 48, 9502–, https://doi.org/10.1029/2011WR011524, 2012.
	"""
	
	def __init__( self , distY = None , distX = None ):##{{{
		"""
		Initialisation of MRec.
		
		Parameters
		----------
		distY: describe the distribution of reference, see QM.
			If None, rv_histogram is used
		distX: describe the distribution of biased dataset, see QM
			If None, rv_histogram is used
		
		"""
		self._qmX0 = None
		self._qmX1 = None
		self._qmY0 = None
		self._S_CY0g  = None
		self._Si_CX0g = None
		self._re_un_mat = None
		self.n_features = 0
		self._distY = distY
		self._distX = distX
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
		if Y0.ndim == 1 : Y0 = Y0.reshape(-1,1)
		if X0.ndim == 1 : X0 = X0.reshape(-1,1)
		if X1.ndim == 1 : X1 = X1.reshape(-1,1)
		self.n_features = Y0.shape[1]
		
		
		## Kind of variables
		if self._distY is None: self._distY = [rv_histogram for _ in range(self.n_features)]
		if self._distX is None: self._distX = [rv_histogram for _ in range(self.n_features)]
		
		## Transform into Gaussian data
		self._qmY0 = QM( distY0 = sc.norm(0,1) , distX0 = self._distY )
		self._qmX0 = QM( distY0 = sc.norm(0,1) , distX0 = self._distX )
		self._qmX1 = QM( distY0 = sc.norm(0,1) , distX0 = self._distX )
		self._qmY0.fit( None , Y0 )
		self._qmX0.fit( None , X0 )
		self._qmX1.fit( None , X1 )
		Y0g = self._qmY0.predict(Y0)
		X0g = self._qmX0.predict(X0)
		X1g = self._qmX1.predict(X1)
		
		## Correlation matrix
		CY0g = np.corrcoef( Y0g.T )
		CX0g = np.corrcoef( X0g.T )
		
		## Squareroot matrix
		a_CY0g,d_CY0g,_ = np.linalg.svd(CY0g)
		self._S_CY0g = a_CY0g @ np.diag(np.sqrt(d_CY0g)) @ a_CY0g.T
		
		a_CX0g,d_CX0g,_ = np.linalg.svd(CX0g)
		self._Si_CX0g = a_CX0g @ np.diag( np.power(d_CX0g,-0.5) ) @ a_CX0g.T
		
		## Decor-recor-relation
		self._re_un_mat = self._S_CY0g @ self._Si_CX0g
		X0_recor = np.transpose( self._re_un_mat @ X0g.T )
		X1_recor = np.transpose( self._re_un_mat @ X1g.T )
		
		## Final QM
		self._qmY0 = QM( distY0 = self._distY , distX0 = sc.norm )
		self._qmY0.fit( Y0 , X0_recor )
	##}}}
	
	def predict( self , X1 , X0 = None ):##{{{
		"""
		Perform the bias correction
		Return Z1 if X0 is None, else return a tuple Z1,Z0
		
		Parameters
		----------
		X1 : np.array[ shape = (n_sample,n_features) ]
			Array of value to be corrected in projection period
		X0 : np.array[ shape = (n_sample,n_features) ] or None
			Array of value to be corrected in calibration period, optional
		
		Returns
		-------
		Z1 : np.array[ shape = (n_sample,n_features) ]
			Return an array of correction in projection period
		Z0 : np.array[ shape = (n_sample,n_features) ] or None
			Return an array of correction in calibration period
		"""
		if X0 is not None and X0.ndim == 1 : X0 = X0.reshape(-1,1)
		if X1.ndim == 1 : X1 = X1.reshape(-1,1)
		
		X1g = self._qmX1.predict(X1)
		X1_recor = np.transpose( self._re_un_mat @ X1g.T )
		Z1 = self._qmY0.predict(X1_recor)
		
		Z0 = None
		if X0 is not None:
			X0g = self._qmX0.predict(X0)
			X0_recor = np.transpose( self._re_un_mat @ X0g.T )
			Z0 = self._qmY0.predict(X0_recor)
			return Z1,Z0
		return Z1
	##}}}


