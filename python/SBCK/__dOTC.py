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
from .__OTC                       import OTC



###########
## Class ##
###########

class dOTC:
	"""
	SBCK.dOTC
	=========
	
	Description
	-----------
	Dynamical Optimal Transport bias Corrector, taking account of an evolution of the distribution. see [1]
	
	References
	----------
	[1] Robin, Y., Vrac, M., Naveau, P., Yiou, P.: Multivariate stochastic bias corrections with optimal transport, Hydrol. Earth Syst. Sci., 23, 773–786, 2019, https://doi.org/10.5194/hess-23-773-2019
	"""
	
	def _eps_cholesky( self , M , nit = 200 ): #{{{
		MC = None
		try:
			MC = np.linalg.cholesky(M)
		except:
			MC = None
		eps = 0
		if MC is None:
			eps = min( 1e-9 , np.abs(np.diagonal(M)).min() )
			if eps >= 0 and eps <= 0:
				eps = 1e-9
			it = 0
			while MC is None and it < nit:
				perturb = np.identity( M.shape[0] ) * eps
				try:
					MC = np.linalg.cholesky( M + perturb )
				except:
					MC = None
				eps = 2 * eps
				nit += 1
		return MC
	#}}}
	
	def __init__( self , bin_width = None , bin_origin = None , cov_factor = "std" , ot = OTNetworkSimplex() ):##{{{
		"""
		Initialisation of dynamical Optimal Transport bias Corrector.
		
		Parameters
		----------
		bin_width  : np.array[ shape = (n_features) ] or None
			Lenght of bins, see SBCK.tools.SparseHist. If None, bin_width is estimated during fit.
		bin_origin : np.array[ shape = (n_features) ] or None
			Corner of one bin, see SBCK.tools.SparseHist. If None, np.repeat( 0 , n_features ) is used.
		cov_factor : str or np.array[ shape = (n_features,n_features) ]
			Correction factor during transfer of the evolution between X0 and X1 to Y0
				"cholesky" => compute the cholesky factor
				"std"      => compute the standard deviation factor
				other str  => identity is used
		ot         : OT*Solver*
			A solver for Optimal transport, default is SBCK.tools.OTNetworkSimplex()
		
		Attributes
		----------
		otc   : SBCK.OTC
			OTC corrector between X1 and the estimation of Y1
		"""
		self._cov_factor_str = cov_factor
		self._cov_factor = None if type(cov_factor) == str else cov_factor
		self.bin_width  = bin_width
		self.bin_origin = bin_origin
		self._otcX0Y0   = None
		self.otc         = None
		self._ot         = ot
	##}}}
	
	def fit( self , Y0 , X0 , X1 ):##{{{
		"""
		Fit the dOTC model to perform non-stationary bias correction during period 1. For period 0, see OTC
		
		Parameters
		----------
		Y0 : np.array[ shape = (n_samples,n_features) ]
			Reference dataset during period 0
		X0 : np.array[ shape = (n_samples,n_features) ]
			Biased dataset during period 0
		X1	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset during period 1
		"""
		## Set the covariance factor correction
		
		if Y0.ndim == 1: Y0 = Y0.reshape(-1,1)
		if X0.ndim == 1: X0 = X0.reshape(-1,1)
		if X1.ndim == 1: X1 = X1.reshape(-1,1)
		
		if self._cov_factor is None:
			if self._cov_factor_str in ["std" , "cholesky"]:
				if Y0.shape[1] == 1:
					try:
						self._cov_factor = np.std( Y0 ) / np.std( X0 )
					except:
						self._cov_factor = 1
				elif self._cov_factor_str == "cholesky":
					fact0 = self._eps_cholesky( np.cov( Y0 , rowvar = False ) )
					fact1 = self._eps_cholesky( np.cov( X0 , rowvar = False ) )
					self._cov_factor = np.dot( fact0 , np.linalg.inv( fact1 ) )
				else:
					fact0 = np.std( Y0 , axis = 0 )
					fact1 = np.std( X0 , axis = 0 )
					self._cov_factor = np.diag( fact0 / fact1 )
			else:
				self._cov_factor = np.identity(Y0.shape[1])
		
		self.bin_width = self.bin_width if self.bin_width is not None else bin_width_estimator( [Y0,X0,X1] )
		
		
		## Optimal plan
		otcY0X0 = OTC( self.bin_width , self.bin_origin , ot = self._ot )
		otcX0X1 = OTC( self.bin_width , self.bin_origin , ot = self._ot )
		otcY0X0.fit( X0 , Y0 )
		otcX0X1.fit( X1 , X0 )
		self._otcX0Y0 = OTC( self.bin_width , self.bin_origin , ot = self._ot )
		self._otcX0Y0.fit(Y0,X0)
		
		## Estimation of Y1
		yX0 = otcY0X0.predict(Y0)
		yX1 = otcX0X1.predict(yX0)
		motion = yX1 - yX0
		motion = np.apply_along_axis( lambda x : np.dot( self._cov_factor , x ) , 1 , motion )
		Y1 = Y0 + motion
		
		## Optimal plan for correction
		self.otc = OTC( self.bin_width , self.bin_origin , ot = self._ot )
		self.otc.fit( Y1 , X1 )
	##}}}
	
	def predict( self , X1 , X0 = None ):##{{{
		"""
		Perform the bias correction
		Return Z1 if X0 is None, else return a tuple Z1,Z0
		
		Note: Only the center of the bins associated to the corrected points are
		returned, but all corrections of the form:
		>> dotc.predict(X1) + np.random.uniform( low = - dotc.bin_width / 2 , high = dotc.bin_width / 2 , size = X1.shape[0] )
		are equivalent for dOTC.
		
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
		Z1 = self.otc.predict( X1 )
		if X0 is not None:
			Z0 = self._otcX0Y0.predict(X0)
			return Z1,Z0
		return Z1
	##}}}
	



