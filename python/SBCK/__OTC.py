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
from .tools.__tools_cpp           import SparseHist
from .tools.__bin_width_estimator import bin_width_estimator
from .tools.__OT                  import OTNetworkSimplex
from .tools.__OT                  import OTSinkhornLogDual


###########
## Class ##
###########

class OTC:
	"""
	SBCK.OTC
	========
	
	Description
	-----------
	Optimal Transport bias Corrector, see [1]
	
	References
	----------
	[1] Robin, Y., Vrac, M., Naveau, P., Yiou, P.: Multivariate stochastic bias corrections with optimal transport, Hydrol. Earth Syst. Sci., 23, 773–786, 2019, https://doi.org/10.5194/hess-23-773-2019
	"""
	
	def __init__( self , bin_width = None , bin_origin = None , ot = OTNetworkSimplex() ):##{{{
		"""
		Initialisation of Optimal Transport bias Corrector.
		
		Parameters
		----------
		bin_width  : np.array( [shape = (n_features) ] )
			Lenght of bins, see SBCK.SparseHist. If is None, it is estimated during the fit
		bin_origin : np.array( [shape = (n_features) ] )
			Corner of one bin, see SBCK.SparseHist. If is None, np.repeat(0,n_features) is used
		ot         : OT*Solver*
			A solver for Optimal transport, default is OTSinkhornLogDual()
		
		Attributes
		----------
		muY	: SBCK.SparseHist
			Multivariate histogram of references
		muX	: SBCK.SparseHist
			Multivariate histogram of biased dataset
		"""
		
		self.muX = None
		self.muY = None
		self.bin_width  = bin_width
		self.bin_origin = bin_origin
		self._plan       = None
		self._ot         = ot
	##}}}
	
	def fit( self , Y0 , X0 ):##{{{
		"""
		Fit the OTC model
		
		Parameters
		----------
		Y0	: np.array[ shape = (n_samples,n_features) ]
			Reference dataset
		X0	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset
		"""
		
		## Sparse Histogram
		self.bin_width  = np.array( [self.bin_width ] ).ravel() if self.bin_width  is not None else bin_width_estimator( [Y0,X0] )
		self.bin_origin = np.array( [self.bin_origin] ).ravel() if self.bin_origin is not None else np.zeros( self.bin_width.size )
		
		self.bin_width  = np.array( [self.bin_width] ).ravel()
		self.bin_origin = np.array( [self.bin_origin] ).ravel()
		
		self.muY = SparseHist( Y0 , bin_width = self.bin_width , bin_origin = self.bin_origin )
		self.muX = SparseHist( X0 , bin_width = self.bin_width , bin_origin = self.bin_origin )
		
		
		## Optimal Transport
		self._ot.fit( self.muX , self.muY )
		if not self._ot.state:
			print( "Warning: Error in network simplex, try SinkhornLogDual" )
			self._ot = OTSinkhornLogDual()
			self._ot.fit( self.muX , self.muY )
		
		## 
		self._plan = np.copy( self._ot.plan() )
		self._plan = ( self._plan.T / self._plan.sum( axis = 1 ) ).T
	##}}}
	
	def predict( self , X0 ):##{{{
		"""
		Perform the bias correction.
		
		Note: Only the center of the bins associated to the corrected points are
		returned, but all corrections of the form:
		>> otc.predict(X0) + np.random.uniform( low = - otc.bin_width / 2 , high = otc.bin_width / 2 , size = X0.shape[0] )
		are equivalent for OTC.
		
		Parameters
		----------
		X0  : np.array[ shape = (n_samples,n_features) ]
			Array of values to be corrected
		
		Returns
		-------
		Z0 : np.array[ shape = (n_samples,n_features) ]
			Return an array of correction
		"""
		indx = self.muX.argwhere(X0)
		indy = np.zeros_like(indx)
		for i,ix in enumerate(indx):
			indy[i] = np.random.choice( range(self.muY.size) , p = self._plan[ix,:] )
		return self.muY.c[indy,:]
	##}}}



