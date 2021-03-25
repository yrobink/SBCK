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
from .__CDFt          import CDFt
from .tools.__shuffle import SchaakeShuffle


###########
## Class ##
###########

class ECBC(CDFt):
	"""
	SBCK.ECBC
	=========
	
	Description
	-----------
	This class implements the method Empirical Copula Bias Correction discribed in [1].
	
	References
	----------
	[1] Vrac, M. and P. Friederichs, 2015: Multivariate—Intervariable, Spatial, and Temporal—Bias Correction. J. Climate, 28, 218–237, https://doi.org/10.1175/JCLI-D-14-00059.1
	"""
	
	def __init__( self , **kwargs ):##{{{
		"""
		Initialisation of ECBC
		
		Parameters
		----------
		
		**kwargs: Any named arguments
			All are passed to SBCK.CDFt
		
		Attributes
		----------
		"""
		CDFt.__init__( self , **kwargs )
		self._ss = SchaakeShuffle()
	##}}}
	
	def fit( self , Y0 , X0 , X1 = None ):##{{{
		"""
		Fit ECBC
		
		Parameters
		----------
		Y0	: np.array[ shape = (n_samples,n_features) ]
			Reference dataset during calibration period
		X0	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset during calibration period
		X1	: np.array[ shape = (n_samples,n_features) ] or None
			Biased dataset during projection period. Can be None to use as a stationary bias correction method
		"""
		CDFt.fit( self , Y0 , X0 , X1 )
		self._ss.fit(Y0)
	##}}}
	
	def predict( self , X1 , X0 = None ):##{{{
		"""
		Perform the bias correction
		
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
			Return an array of correction in calibration period, or None
		"""
		
		Z = CDFt.predict( self , X1 , X0 )
		if X0 is not None:
			Z1,Z0 = Z
			Z1 = self._ss.predict(Z1)
			Z0 = self._ss.predict(Z0)
			return Z1,Z0
		Z1 = Z
		Z1 = self._ss.predict(Z1)
		return Z1
	##}}}



