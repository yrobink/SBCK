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
from .__OTC import OTC
from .tools.__Shift import Shift


###########
## Class ##
###########

class TSMBC:
	"""
	SBCK.TSMBC
	==========
	
	Description
	-----------
	Time Shifted Multivariate Bias Correction.
	
	References
	----------
	[1] Robin, Y. and Vrac, M.: Is time a variable like the others in multivariate statistical downscaling and bias correction?, Earth Syst. Dynam. Discuss. [preprint], https://doi.org/10.5194/esd-2021-12, in review, 2021.
	"""
	
	def __init__( self , lag , bc_method = OTC , method = "row" , ref = "middle" , **kwargs ):##{{{
		"""
		Initialisation of TSMBC.
		
		Parameters
		----------
		lag       : integer
			Time lag of the shift
		bc_method : An class of SBCK
			bias correction method used, default is SBCK.OTC
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
	
	## Properties {{{
	
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
	
	def fit( self , Y0 , X0 ):##{{{
		"""
		Fit of the bc_method model on shifted X0 and Y0
		
		Parameters
		----------
		Y0	: np.array[ shape = (n_samples,n_features) ]
			Reference dataset
		X0	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset
		"""
		Xs = self.shift.transform(X0)
		Ys = self.shift.transform(Y0)
		self.bc_method.fit( Ys , Xs )
	##}}}
	
	def predict( self , X0 ):##{{{
		"""
		Perform the bias correction
		
		Parameters
		----------
		X0  : np.array[ shape = (n_samples,n_features) ]
			Array of values to be corrected
		
		Returns
		-------
		Z0 : np.array[ shape = (n_samples,n_features) ]
			Return an array of correction
		"""
		Xs = self.shift.transform(X0)
		return self.shift.inverse( self.bc_method.predict(Xs) )
	##}}}


