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
from .__QM import QM


###########
## Class ##
###########

class QDM:
	"""
	SBCK.QDM
	========
	
	Description
	-----------
	QDM Bias correction method, see [1]
	
	References
	----------
	[1] Cannon, A. J., Sobie, S. R., and Murdock, T. Q.: Bias correction of simulated precipitation by quantile mapping: how well do methods preserve relative changes in quantiles and extremes?, J. Climate, 28, 6938–6959, https://doi.org/10.1175/JCLI-D-14- 00754.1, 2015.
	"""
	
	def __init__( self , delta = "additive" , **kwargs ):##{{{
		"""
		Initialisation of QDM.
		
		Parameters
		----------
		delta : str or tuple
			Delta method: "additive" or "multiplicative". It is possible to pass custom delta function with a tuple where the first element is the transform ( e.g. np.add or np.multiply) and the second one its inverse (e.g. np.subtract or np.divide)
		kwargs: any named arguments
			QDM call the QM method, any arguments of kwargs are passed to QM.
		
		"""
		self._delta_method  = np.add
		self._idelta_method = np.subtract
		if delta == "multiplicative":
			self._delta_method  = np.multiply
			self._idelta_method = np.divide
		if isinstance(delta,(list,tuple)):
			self._delta_method  = delta[0]
			self._idelta_method = delta[1]
		
		self._delta  = None
		self._qmX1Y0 = None
		self._qm_args = kwargs
	##}}}
	
	def fit( self , Y0 , X0 , X1 ):##{{{
		"""
		Fit the MBCn model
		
		Parameters
		----------
		Y0 : np.array[ shape = (n_samples,n_features) ]
			Reference dataset during calibration period
		X0 : np.array[ shape = (n_samples,n_features) ]
			Biased dataset during calibration period
		X1	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset during projection period
		"""
		self._qmX0Y0 = QM(**self._qm_args)
		self._qmX0Y0.fit(Y0,X0)
		qmX1X0 = QM(**self._qm_args)
		qmX1X0.fit( X0 , X1 )
		self._delta = self._idelta_method( X1 , qmX1X0.predict(X1) )
		self._qmX1Y0 = QM(**self._qm_args)
		self._qmX1Y0.fit( Y0 , X1 )
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
		Z0 = None if X0 is None else self._qmX0Y0.predict(X0)
		Z1 = self._delta_method( self._qmX1Y0.predict(X1) , self._delta )
		if Z0 is not None:
			return Z1,Z0
		return Z1
	##}}}


