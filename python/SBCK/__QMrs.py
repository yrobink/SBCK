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
from .__QM            import QM
from .tools.__shuffle import SchaakeShuffleRef


###########
## Class ##
###########

class QMrs(QM):
	"""
	SBCK.QMrs
	=========
	
	Description
	-----------
	Quantile Mapping bias corrector with multivariate rankshuffle, as described in [1]
	
	References
	----------
	[1] Vrac, M.: Multivariate bias adjustment of high-dimensional climate simulations: the Rank Resampling for Distributions and Dependences (R2 D2 ) bias correction, Hydrol. Earth Syst. Sci., 22, 3175–3196, https://doi.org/10.5194/hess-22-3175-2018, 2018.
	"""
	
	
	def __init__( self , refs = [0] , **kwargs ):##{{{
		"""
		Initialisation of Quantile Mapping bias corrector.
		
		Parameters
		----------
		refs     : list
			Index of reference for rankshuffle, see SBCK.tools.SchaakeShuffleRef
		**kwargs : see SBCK.QM
			All others arguments are passed to SBCK.QM class.
		"""
		QM.__init__( self , **kwargs )
		self._refs = refs
		self._ssr  = SchaakeShuffleRef( refs[0] )
	##}}}
	
	def fit( self , Y0 , X0 ):##{{{
		"""
		Fit of the quantile mapping model
		
		Parameters
		----------
		Y0	: np.array[ shape = (n_samples,n_features) ]
			Reference dataset
		X0	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset
		"""
		QM.fit( self , Y0 , X0 )
		self._ssr.fit(Y0)
	##}}}
	
	def predict( self , X0  ):##{{{
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
		Zu = QM.predict( self , X0 )
		Z0 = np.zeros( Zu.shape + (len(self._refs),) )
		for i,r in enumerate(self._refs):
			self._ssr._ref = r
			Z0[:,:,i] = self._ssr.predict(X0)
		return Z0.squeeze()
	##}}}

