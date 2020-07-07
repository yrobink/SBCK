# -*- coding: utf-8 -*-

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the SBCK (Statistical    ##
## Bias Correction Kit). This library makes it possible to perform bias         ##
## correction with non parametric methods, and give some metrics between Sparse ##
## Histogram is high dimensions.                                                ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## SBCK (Statistical Bias Correction Kit). Cette librairie permet d'appliquer   ##
## une correction de biais avec des méthodes non paramétriques, et propose      ##
## diverses metrique entre Histograme Sparse en haute dimension.                ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Original author  : Mathieu Vrac                                              ##
## Contact          : mathieu.vrac@lsce.ipsl.fr                                 ##
##                                                                              ##
## Notes   : CDFt is the re-implementation of the function CDFt of R package    ##
##           "CDFt" developped by Mathieu Vrac, available at                    ##
##           https://cran.r-project.org/web/packages/CDFt/index.html            ##
##           This code is governed by the CeCILL-C license with the             ##
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
			Numerical tolerance, default 1e-3
		bin_width : np.array[ shape = (n_features) ]
			Lenght of bins for each margins. If None, length of bins are estimating during fit.
		
		"""
		self.n_features = kwargs.get("n_features")
		self._tol = kwargs.get("tol") if kwargs.get("tol") is not None else 1e-3
		
		self._distY0 = _Dist( dist = kwargs.get("distY0") , kwargs = kwargs.get("kwargsY0") )
		self._distY1 = _Dist( dist = kwargs.get("distY1") , kwargs = kwargs.get("kwargsY1") )
		self._distX0 = _Dist( dist = kwargs.get("distX0") , kwargs = kwargs.get("kwargsX0") )
		self._distX1 = _Dist( dist = kwargs.get("distX1") , kwargs = kwargs.get("kwargsX1") )
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
					Y1 = self._infer_Y1( Y0uni , X0uni , X1uni , i )
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

