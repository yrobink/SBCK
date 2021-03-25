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

from .tools.__Dist import _Dist

###########
## Class ##
###########

class QM:
	"""
	SBCK.QM
	=======
	
	Description
	-----------
	Quantile Mapping bias corrector, see e.g. [1,2,3]. The implementation proposed here is generic, and can use
	scipy.stats to fit a parametric distribution, or can use a frozen distribution.
	
	Example
	-------
	```
	## Start with a reference / biased dataset, noted Y,X, from normal distribution:
	X = np.random.normal( loc = 0 , scale = 2   , size = 1000 )
	Y = np.random.normal( loc = 5 , scale = 0.5 , size = 1000 )
	
	## Generally, we do not know the distribution of X and Y, and we use the empirical quantile mapping:
	qm_empiric = QM( distY0 = SBCK.tools.rv_histogram , distX0 = SBCK.tools.rv_histogram ) ## = QM(), default
	qm_empiric.fit(Y,X)
	Z_empiric = qm_empiric.predict(X) ## Z is the correction in a non parametric way
	
	## But we can know that X and Y follow a Normal distribution, without knowing the parameters:
	qm_normal = QM( distY0 = scipy.stats.norm , distX0 = scipy.stats.norm )
	qm_normal.fit(Y,X)
	Z_normal = qm_normal.predict(X)
	
	## And finally, we can know the law of Y, and it is usefull to freeze the distribution:
	qm_freeze = QM( distY0 = scipy.stats.norm( loc = 5 , scale = 0.5 ) , distX0 = scipy.stats.norm )
	qm_freeze.fit(Y,X) ## = qm_freeze.fit(None,X) because Y is not used
	Z_freeze = qm_freeze.predict(X)
	```
	
	References
	----------
	[1] Panofsky, H. A. and Brier, G. W.: Some applications of statistics to meteorology, Mineral Industries Extension Services, College of Mineral Industries, Pennsylvania State University, 103 pp., 1958.
	[2] Wood, A. W., Leung, L. R., Sridhar, V., and Lettenmaier, D. P.: Hydrologic Implications of Dynamical and Statistical Approaches to Downscaling Climate Model Outputs, Clim. Change, 62, 189–216, https://doi.org/10.1023/B:CLIM.0000013685.99609.9e, 2004.
	[3] Déqué, M.: Frequency of precipitation and temperature extremes over France in an anthropogenic scenario: Model results and statistical correction according to observed values, Global Planet. Change, 57, 16–26, https://doi.org/10.1016/j.gloplacha.2006.11.030, 2007.
	"""
	
	def __init__( self , **kwargs ):##{{{
		"""
		Initialisation of Quantile Mapping bias corrector. All arguments must be named.
		
		Parameters
		----------
		distY0 : A statistical distribution from scipy.stats or SBCK.tools.rv_*
			The distribution of references.
		distX0 : A statistical distribution from scipy.stats or SBCK.tools.rv_*
			The distribution of biased dataset.
		kwargsY0 : dict
			Arguments passed to distY0
		kwargsX0 : dict
			Arguments passed to distX0
		n_features: None or integer
			Numbers of features, optional because it is determined during fit if X0 and Y0 are not None.
		tol : float
			Numerical tolerance, default 1e-3
		"""
		self.n_features = kwargs.get("n_features")
		self._tol = kwargs.get("tol") if kwargs.get("tol") is not None else 1e-3
		
		self._distY0 = _Dist( dist = kwargs.get("distY0") , kwargs = kwargs.get("kwargsY0") )
		self._distX0 = _Dist( dist = kwargs.get("distX0") , kwargs = kwargs.get("kwargsX0") )
	##}}}
	
	def fit( self , Y0 , X0 ):##{{{
		"""
		Fit the QM model
		
		Parameters
		----------
		Y0	: np.array[ shape = (n_samples,n_features) ]
			Reference dataset
		X0	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset
		"""
		## Reshape data in form [n_samples,n_features]
		if Y0 is not None and Y0.ndim == 1 : Y0 = Y0.reshape(-1,1)
		if X0 is not None and X0.ndim == 1 : X0 = X0.reshape(-1,1)
		if self.n_features is None:
			if Y0 is None and X0 is None:
				print( "n_features must be set during initialization if Y0 = X0 = None" )
			elif Y0 is not None: self.n_features = Y0.shape[1]
			else: self.n_features = X0.shape[1]
		
		## 
		self._distY0.set_features(self.n_features)
		self._distX0.set_features(self.n_features)
		
		## Fit
		for i in range(self.n_features):
			if Y0 is not None: self._distY0.fit( Y0[:,i] , i )
			else : self._distY0.fit( None , i )
			if X0 is not None: self._distX0.fit( X0[:,i] , i )
			else : self._distX0.fit( None , i )
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
		if X0.ndim == 1 : X0 = X0.reshape(-1,1)
		Z0 = np.zeros_like(X0)
		for i in range(self.n_features):
			cdf = self._distX0.law[i].cdf(X0[:,i])
			cdf[np.logical_not(cdf < 1)] = 1 - self._tol
			cdf[np.logical_not(cdf > 0)] = self._tol
			Z0[:,i] = self._distY0.law[i].ppf( cdf )
		
		return Z0
	##}}}


