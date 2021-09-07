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
from .__dOTC import dOTC
from .tools.__Shift import Shift


###########
## Class ##
###########

class dTSMBC:
	"""
	SBCK.dTSMBC
	===========
	
	Description
	-----------
	Time Shifted Multivariate Bias Correction where observations are unknown.
	
	References
	----------
	[1] Robin, Y. and Vrac, M.: Is time a variable like the others in multivariate statistical downscaling and bias correction?, Earth Syst. Dynam. Discuss. [preprint], https://doi.org/10.5194/esd-2021-12, in review, 2021.
	"""
	
	def __init__( self , lag , bc_method = dOTC , method = "row" , ref = "middle" , **kwargs ):##{{{
		"""
		Initialisation of dTSMBC.
		
		Parameters
		----------
		lag       : integer
			Time lag of the shift
		bc_method : An element of SBCK
			bias correction method used, default is SBCK.dOTC()
		method    : string
			inverse method for shift, see SBCK.tools.Shift
		ref       : integer
			Reference columns/rows for inverse, see SBCK.tools.Shift, default is 0.5 * (lag+1)
		**kwargs  : arguments of bc_method
		
		Attributes
		----------
		bc_method : An element of SBCK
			Bias correction method
		shift     : Shift class
			class used to shift and un-shift data
		"""
		self.bc_method = bc_method(**kwargs)
		if ref == "middle": ref = int(0.5*(lag+1))
		self.shift     = Shift( lag , method , ref )
	##}}}
	
	## Methods and properties ##{{{
	
	@property
	def ref(self):
		return self.shift.ref
	
	@ref.setter
	def ref( self , _ref ):
		self.shift.ref = _ref
	
	@property
	def method(self):
		return self.shift.method
	
	@method.setter
	def method( self , _method ):
		self.shift.method = _method
	##}}}
	
	def fit( self , Y0 , X0 , X1 ):##{{{
		"""
		Fit of the bc_method model on shifted X1, with learning shifted pair of Y0 and X0
		
		Parameters
		----------
		Y0	: np.array[ shape = (n_samples,n_features) ]
			Reference dataset on learning part
		X0	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset on learning part
		X1	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset on projection part
		"""
		Y0s = self.shift.transform(Y0)
		X0s = self.shift.transform(X0)
		X1s = self.shift.transform(X1)
		self.bc_method.fit( Y0s , X0s , X1s )
	##}}}
	
	def predict( self , X1 , X0 = None ):##{{{
		"""
		Perform the bias correction of the shifted X1, and return the unshift correction
		Return Z1 if X0 is None, else return a tuple Z1,Z0
		
		Parameters
		----------
		X1  : np.array[ shape = (n_samples,n_features) ]
			Array of values to be corrected in projection period
		X0  : np.array[ shape = (n_samples,n_features) ] or None
			Array of values to be corrected in calibration period
		
		Returns
		-------
		Z1 : np.array[ shape = (n_sample,n_features) ]
			Return an array of correction in projection period
		Z0 : np.array[ shape = (n_sample,n_features) ] or None
			Return an array of correction in calibration period, or None
		"""
		X1s = self.shift.transform(X1)
		if X0 is not None:
			X0s = self.shift.transform(X0)
			Z1s,Z0s = self.bc_method.predict(X1s,X0s)
			Z1 = self.shift.inverse(Z1s)
			Z0 = self.shift.inverse(Z0s)
			return Z1,Z0
		
		Z1s = self.bc_method.predict(X1s)
		Z1 = self.shift.inverse(Z1s)
		return Z1
	##}}}

