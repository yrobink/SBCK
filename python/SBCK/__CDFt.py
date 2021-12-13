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

##################################################################################
##################################################################################
##                                                                              ##
## Original author  : Mathieu Vrac                                              ##
## Contact          : mathieu.vrac@lsce.ipsl.fr                                 ##
##                                                                              ##
## Notes   : CDFt is the re-implementation of the function CDFt of R package    ##
##           "CDFt" developped by Mathieu Vrac, available at                    ##
##           https://cran.r-project.org/web/packages/CDFt/index.html            ##
##           This code is governed by the GNU-GPL3 license with the             ##
##           authorization of Mathieu Vrac                                      ##
##                                                                              ##
##################################################################################
##################################################################################


###############
## Libraries ##
###############

import numpy       as np
import scipy.stats as sc
import scipy.interpolate as sci

from .tools.__Dist import _Dist


###########
## Class ##
###########

class CDFt:
	"""
	SBCK.CDFt
	=========
	
	Description
	-----------
	Quantile Mapping bias corrector, taking account of an evolution of the distribution, see [1].
	
	References
	----------
	[1] Michelangeli, P.-A., Vrac, M., and Loukos, H.: Probabilistic downscaling approaches: Application to wind cumulative distribution functions, Geophys. Res. Lett., 36, L11708, https://doi.org/10.1029/2009GL038401, 2009.
	
	Notes
	-----
	CDFt is the re-implementation of the function CDFt of R package "CDFt" developped by Mathieu Vrac, available at
	https://cran.r-project.org/web/packages/CDFt/index.htmm
	"""
	def __init__( self , **kwargs ):##{{{
		"""
		Initialisation of CDFt bias corrector. All arguments must be named.
		
		Parameters
		----------
		distY0 : A statistical distribution from scipy.stats or SBCK.tools.rv_*
			The distribution of references in calibration period. Default is rv_histogram.
		distX0 : A statistical distribution from scipy.stats or SBCK.tools.rv_*
			The distribution of biased dataset in calibration period. Default is rv_histogram.
		distY1 : A statistical distribution from scipy.stats or SBCK.tools.rv_*
			The distribution of references in projection period. Default is rv_histogram, and Y1 is estimated during fit
		distX1 : A statistical distribution from scipy.stats or SBCK.tools.rv_*
			The distribution of biased dataset in projection period. Default is rv_histogram.
		kwargsY0 : dict
			Arguments passed to distY0
		kwargsX0 : dict
			Arguments passed to distX0
		kwargsY1 : dict
			Arguments passed to distY1
		kwargsX1 : dict
			Arguments passed to distX1
		n_features: None or integer
			Numbers of features, optional because it is determined during fit if X0 and Y0 are not None.
		tol : float
			Numerical tolerance, default 1e-6
		
		"""
		self.n_features = kwargs.get("n_features")
		self._tol   = kwargs.get("tol") if kwargs.get("tol") is not None else 1e-6
		self._dsupp = kwargs.get("dsupp") if kwargs.get("dsupp") is not None else 1000
		self._samples_Y1 = kwargs.get("samples_Y1") if kwargs.get("samples_Y1") is not None else 10000
		
		self._distY0 = _Dist( dist = kwargs.get("distY0") , kwargs = kwargs.get("kwargsY0") )
		self._distY1 = _Dist( dist = kwargs.get("distY1") , kwargs = kwargs.get("kwargsY1") )
		self._distX0 = _Dist( dist = kwargs.get("distX0") , kwargs = kwargs.get("kwargsX0") )
		self._distX1 = _Dist( dist = kwargs.get("distX1") , kwargs = kwargs.get("kwargsX1") )
		self._scale_tail = kwargs.get("scale_tail")
	##}}}
	
	def fit( self , Y0 , X0 , X1 ):##{{{
		"""
		Fit of CDFt model
		
		Parameters
		----------
		Y0	: np.array[ shape = (n_samples,n_features) ]
			Reference dataset during calibration period
		X0	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset during calibration period
		X1	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset during projection period
		
		Note
		----
		The fit is performed margins by margins (without taking into account the dependance structure, see R2D2 or dOTC)
		"""
		
		## Reshape data in matrix form
		if Y0 is not None and Y0.ndim == 1 : Y0 = Y0.reshape(-1,1)
		if X0 is not None and X0.ndim == 1 : X0 = X0.reshape(-1,1)
		if X1 is not None and X1.ndim == 1 : X1 = X1.reshape(-1,1)
		
		## Find n_features
		if self.n_features is None:
			if Y0 is None and X0 is None and X1 is None:
				print( "n_features must be set during initialization if Y0 = X0 = X1 = None" )
			elif Y0 is not None: self.n_features = Y0.shape[1]
			elif X0 is not None: self.n_features = X0.shape[1]
			else:                self.n_features = X1.shape[1]
		
		## Find laws
		self._distY0.set_features(self.n_features)
		self._distY1.set_features(self.n_features)
		self._distX0.set_features(self.n_features)
		self._distX1.set_features(self.n_features)
		
		## Start fit itself
		for i in range(self.n_features):
			self._distY0.fit( Y0[:,i] , i )
			self._distX0.fit( X0[:,i] , i )
			self._distX1.fit( X1[:,i] , i )
			## Fit Y1
			if self._distY1.is_frozen(i):
				self._distY1.law.append(self._dist.distY1[i])
			else:
				if self._distY0.is_parametric(i) and self._distX0.is_parametric(i) and self._distX1.is_parametric(i):
					Y1 = self._distX1.law[i].ppf( self._distX0.law[i].cdf( self._distY0.law[i].ppf( self._distX1.law[i].cdf(X1[:,i].squeeze()) ) ) )
				else:
					Y0uni = Y0[:,i] if Y0 is not None else self._distY0.law[-1].rvs(10000)
					X0uni = X0[:,i] if X0 is not None else self._distX0.law[-1].rvs(10000)
					X1uni = X1[:,i] if X1 is not None else self._distX1.law[-1].rvs(10000)
					Y1 = self._infer_Y1_new( Y0uni , X0uni , X1uni , i )
				self._distY1.fit( Y1 , i )
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
		
		Note
		----
		The correction is performed margins by margins (without taking into account the dependance structure, see R2D2 or dOTC)
		"""
		if X1.ndim == 1 : X1 = X1.reshape(-1,1)
		Z1 = np.zeros_like(X1)
		for i in range(self.n_features):
			cdf = self._distX1.law[i].cdf(X1[:,i])
			cdf[np.logical_not(cdf < 1)] = 1 - self._tol
			cdf[np.logical_not(cdf > 0)] = self._tol
			Z1[:,i] = self._distY1.law[i].ppf( cdf )
		
		if X0 is not None:
			if X0.ndim == 1 : X0 = X0.reshape(-1,1)
			Z0 = np.zeros_like(X0)
			for i in range(self.n_features):
				cdf = self._distX0.law[i].cdf(X0[:,i])
				cdf[np.logical_not(cdf < 1)] = 1 - self._tol
				cdf[np.logical_not(cdf > 0)] = self._tol
				Z0[:,i] = self._distY0.law[i].ppf( cdf )
			return Z1,Z0
		return Z1
	##}}}
	
	def _infer_Y1( self , Y0 , X0 , X1 , idx ):##{{{
		mY0 = np.mean(Y0)
		mX0 = np.mean(X0)
		
		X0s = X0 + mY0 - mX0
		X1s = X1 + mY0 - mX0
		
		rvY0 = self._distY0.law[idx]
		rvX0 = self._distX0.dist[idx]( *self._distX0.dist[idx].fit( X0s.squeeze()) , **self._distX0.kwargs )
		rvX1 = self._distX1.dist[idx]( *self._distX1.dist[idx].fit( X1s.squeeze()) , **self._distX1.kwargs )
		
		xdiff = abs( np.mean(X1) - np.mean(X0) )
		xdev = 2
		dev_ok = False
		while not dev_ok:
			dev_ok = True
			xmin = min( [ T.min() for T in [Y0,X0,X1] ] ) - xdev * xdiff
			xmax = max( [ T.max() for T in [Y0,X0,X1] ] ) + xdev * xdiff
			
			x = np.linspace( xmin , xmax , 200 )
			cdfY0 = rvY0.cdf(x)
			cdfX0 = rvX0.cdf(x)
			cdfX1 = rvX1.cdf(x)
			cdfY1 = rvY0.cdf(rvX0.ppf(cdfX1))
			
			## Correction of left part
			if Y0.min() < X1s.min():
				i = np.argwhere( x < rvY0.ppf(cdfY1[0]) ).max()
				j = np.argwhere( x < X1s.min() ).max()
				if i < j:
					cdfY1[:(j-i)] = 0
					cdfY1[(j-i):(j+1)] = cdfY0[:(i+1)]
				else:
					cdfY1[:(j+1)] = cdfY0[(i-j):(i+1)]
			
			## Correction of right part
			if cdfY1[-1] < 1:
				i = np.argwhere(x < rvY0.ppf(cdfY1[-1])).max()
				try:
					j = np.argwhere(cdfY1[:-1] == cdfY1[-1]).min()
					dif = min( [x.size - j,x.size - i] )
					cdfY1[j:(j+dif)] = cdfY0[i:(i+dif)]
					if j + dif < x.size: cdfY1[(j+dif):] = 1.
				except:
					xdev *= 2
					dev_ok = False
		
		
		## Find mean deviation
		if False:
			rvY0 = self._lawY0[idx]
			rvX0 = self._lawX0[idx]
			rvX1 = self._lawX1[idx]
			cdfX0diag = rvX0.cdf(rvY0.ppf( rvX1.cdf(x) ))
			test_diag = np.logical_and( cdfX0diag < 1 , cdfX0diag > 0 )
			if np.any( test_diag ):
				print("here")
				Y1_hard = rvX1.ppf(cdfX0diag)
				mean_dev = np.mean(Y1_hard[test_diag] - x[test_diag])
				x -= mean_dev
		
		## Final estimation
		cdfY1_fct = sci.interp1d( cdfY1 , x )
		Y1 = cdfY1_fct( np.random.uniform( size = 10000 , low = cdfY1.min() , high = cdfY1.max() ) )
		
		return Y1
	##}}}
	
	def _infer_Y1_new( self , Y0 , X0 , X1 , idx ):##{{{
		samples_Y1 = self._samples_Y1
		dsupp = self._dsupp
		tol   = self._tol
		
		## Normalization
		mY0 = np.mean(Y0)
		mX0 = np.mean(X0)
		mX1 = np.mean(X1)
		sY0 = np.std(Y0)
		sX0 = np.std(X0)
		
		X0s = (X0 - mX0) * sY0 / sX0 + mY0
		X1s = (X1 - mX1) * sY0 / sX0 + mX1 + mY0 - mX0
		
		## CDF
		rvY0  = self._distY0.law[idx]
		rvX0s = self._distX0.dist[idx]( *self._distX0.dist[idx].fit( X0s.squeeze()) , **self._distX0.kwargs )
		rvX1s = self._distX1.dist[idx]( *self._distX1.dist[idx].fit( X1s.squeeze()) , **self._distX1.kwargs )
		
		## Support
		## Here the support is such that the CDF of Y0, X0s and X1s start from 0
		## and go to 1
		x_min = min([T.min() for T in [Y0,X0s,X1s,X0,X1]])
		x_max = max([T.max() for T in [Y0,X0s,X1s,X0,X1]])
		x_eps = 0.05 * (x_max - x_min)
		x_fac = 1
		x = np.linspace( x_min - x_fac * x_eps , x_max + x_fac * x_eps , dsupp )
		
		def support_test( rv , x ):
			if not abs(rv.cdf(x[0])) < tol:
				return False
			if not abs(rv.cdf(x[-1])-1) < tol:
				return False
			return True
		
		while (not support_test(rvY0,x)) or (not support_test(rvX0s,x)) or (not support_test(rvX1s,x)):
			x_fac *= 2
			x = np.linspace( x_min - x_fac * x_eps , x_max + x_fac * x_eps , dsupp )
		x_fac /= 2
		dsupp = int(dsupp/1.2)
		
		## Loop to check the support
		extend_support = True
		while extend_support:
			
			## Support
			extend_support = False
			dsupp = int(dsupp*1.2)
			x_fac *= 2
			x = np.linspace( x_min - x_fac * x_eps , x_max + x_fac * x_eps , dsupp )
			
			## Inference of the CDF of Y1
			cdfY1 = rvY0.cdf(rvX0s.ppf(rvX1s.cdf(x)))
			
			## Correction of the CDF, we want that the CDF of Y1 start from 0 and goto 1
			if (not abs(cdfY1[0]) < tol):
				## CDF not start at 0
				idx  = np.max(np.argwhere(np.abs(cdfY1[0] - cdfY1) < tol))
				if idx == 0:
					extend_support = True
				else:
					supp_l_X0s = rvX0s.ppf(cdfY1[0]) - rvX0s.ppf(0)
					supp_l_X1s = rvX1s.ppf(cdfY1[0]) - rvX1s.ppf(0)
					supp_l_Y0  = rvY0.ppf(cdfY1[0])  - rvY0.ppf(0)
					supp_l_Y1  = supp_l_Y0 * supp_l_X1s / supp_l_X0s
					if self._scale_tail is not None:
						supp_l_Y1 = self._scale_tail
					if x[idx] - supp_l_Y1 < x[0]:
						extend_support = True
					else:
						idxl = np.argmin(np.abs(x - (x[idx] - supp_l_Y1)))
						cdfY1[:idxl] = 0
						cdfY1[idxl:idx] = rvY0.cdf( np.linspace( rvY0.ppf(0) , rvY0.ppf(cdfY1[idx]) , idx - idxl ) )
			
			if (not abs(cdfY1[-1] - 1) < tol):
				## CDF not finished at 1
				idx = np.min(np.argwhere(np.abs(cdfY1[-1] - cdfY1) < tol))
				if idx == dsupp -1:
					extend_support = True
				else:
					supp_r_X0s = rvX0s.ppf(1) - rvX0s.ppf(cdfY1[-1]) 
					supp_r_X1s = rvX1s.ppf(1) - rvX1s.ppf(cdfY1[-1]) 
					supp_r_Y0  = rvY0.ppf(1)  - rvY0.ppf(cdfY1[-1])  
					supp_r_Y1  = supp_r_Y0 * supp_r_X1s / supp_r_X0s
					if self._scale_tail is not None:
						supp_l_Y1 = self._scale_tail
					if x[idx] + supp_r_Y1 > x[-1]:
						extend_support = True
					else:
						idxr = np.argmin(np.abs(x - (x[idx] + supp_r_Y1)))
						cdfY1[idxr:] = 1
						cdfY1[idx:idxr] = rvY0.cdf( np.linspace( rvY0.ppf(cdfY1[idx]) , rvY0.ppf(1) , idxr - idx ) )
		
		## Draw Y1
		icdfY1 = sci.interp1d( cdfY1 , x )
		hY1    = icdfY1( np.random.uniform( size = samples_Y1 , low = 0 , high = 1 ) )
		
		return hY1
	##}}}
	
