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

from .tools.__bin_width_estimator import bin_width_estimator


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
	def __init__( self , bin_width = None ):##{{{
		"""
		Initialisation of CDFt bias corrector.
		
		Parameters
		----------
		bin_width : np.array[ shape = (n_features) ]
			Lenght of bins for each margins. If None, length of bins are estimating during fit.
		
		"""
		
		self._n_features = 1
		self._diff       = None
		self._dev        = None
		self.bin_width  = bin_width
		self.bins       = None
		self.rvX1        = None
		self.rvY1       = None
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
		if Y0.ndim == 1: Y0 = Y0.reshape( (Y0.size,1) )
		if X0.ndim == 1: X0 = X0.reshape( (X0.size,1) )
		if X1.ndim == 1: X1 = X1.reshape( (X1.size,1) )
		
		self._n_features = Y0.shape[1]
		self._diff = np.zeros( (self._n_features) )
		self._dev = np.repeat( 1. , self._n_features )
		self.bin_width = self.bin_width if self.bin_width is not None else bin_width_estimator( [Y0,X0,X1] )
		self.bins  = [None for i in range(self._n_features)]
		self.rvX1  = [None for i in range(self._n_features)]
		self.rvY1 = [None for i in range(self._n_features)]
		for idx in range(self._n_features):
			self._fit( Y0[:,idx] , X0[:,idx] , X1[:,idx] , idx )
	##}}}
	
	def predict( self , X1 ):##{{{
		"""
		Perform the bias correction
		
		Parameters
		----------
		X1  : np.array[ shape = (n_samples,n_features) ]
			Value or array of value to be corrected
		
		Returns
		-------
		uX1 : np.array[ shape = (n_samples,n_features) ]
			Return an array or a value of correction of X1
		
		Note
		----
		The correction is performed margins by margins (without taking into account the dependance structure, see R2D2 or dOTC)
		"""
		if len(X1.shape) == 1:
			X1 = X1.reshape( (X1.size,1) )
		Z = np.zeros_like(X1)
		for i in range(self._n_features):
			cdfX1 = self.rvX1[i].cdf( X1[:,i] + self._diff[i] )
			cdfX1[ cdfX1 < 0 ] = 0
			cdfX1[ cdfX1 > 1 ] = 1
			Z[:,i] = self.rvY1[i].ppf(cdfX1)
		
		return Z
	##}}}
	
	def _fit( self , Y0 , X0 , X1 , idx ):##{{{
		## center data on mean(Y0)
		self._diff[idx] = np.mean(Y0) - np.mean(X0) 
		X0c = X0 + self._diff[idx] 
		X1c = X1 + self._diff[idx] 
		
		## Bin
		m = np.abs( np.mean( X1 ) - np.mean( X0 ) )
		s = np.std( X1 ) / np.std( X0 )
		s = 1 if s < 1 else s
		
		## Global loop
		dev_ok = False
		while not dev_ok:
			Min = min( [ Z.min() for Z in [Y0,X0,X1] ] ) - m * s * self._dev[idx]
			Max = max( [ Z.max() for Z in [Y0,X0,X1] ] ) + m * s * self._dev[idx]
			self.bins[idx] = np.arange( Min , Max , self.bin_width[idx] )
			
			## Random wariables
			rvY0           = sc.rv_histogram( np.histogram( Y0  , self.bins[idx] ) )
			rvX0           = sc.rv_histogram( np.histogram( X0c , self.bins[idx] ) )
			self.rvX1[idx] = sc.rv_histogram( np.histogram( X1c , self.bins[idx] ) )
			
			## Corrector
			cdfX1 = self.rvX1[idx].cdf( self.bins[idx] )
			cdfY0 = rvY0.cdf( self.bins[idx] )
			cdfX1[ cdfX1 < 0 ] = 0
			cdfX1[ cdfX1 > 1 ] = 1
			cdfY0[ cdfY0 < 0 ] = 0
			cdfY0[ cdfY0 > 1 ] = 1
			cdfY1E = rvY0.cdf( rvX0.ppf( cdfX1 ) )
			cdfY1E[ cdfY1E < 0 ] = 0
			cdfY1E[ cdfY1E > 1 ] = 1
			
			## Correction of left part of CDF
			if Y0.min() < X1c.min():
				Y0q = np.percentile( Y0 , q = 100 * cdfY1E[0] )
				i = np.argwhere( self.bins[idx] <= Y0q ).ravel().max()
				j = np.argwhere( self.bins[idx] < X1c.min() ).ravel().max()
				if i < j:
					cdfY1E[(j-i):j] = cdfY0[:i]
					cdfY1E[:(j-i)] = 0
				else:
					cdfY1E[:j] = cdfY0[(i-j):i]
			
			## Correction of right part of CDF
			kk = np.argwhere( cdfY1E == 1. ).ravel()
			k = cdfY1E.size - 1
			if kk.size > 0:
				k = kk.min()
			if np.abs( cdfY1E[k] - cdfY1E[k-1] ) > 1. / len(self.bins[idx]) or k == cdfY1E.size - 1:
				Y0q = np.percentile( Y0 , q = 100 * cdfY1E[k-1] )
				i = np.argwhere( Y0q <= self.bins[idx] ).ravel().min()
				j = np.argwhere( cdfY1E == cdfY1E[k-1] ).ravel().min()
				if j > 0:
					diff = self.bins[idx].size - max( j , i )
					cdfY1E[j:(j+diff)] = cdfY0[i:(i+diff)]
					if j + diff < self.bins[idx].size:
						cdfY1E[(j+diff):] = 1
					dev_ok = True
				else:
					self._dev[idx] *= 2
			else:
				dev_ok = True
		
		## Random variable of the estimation of Y1
		self.rvY1[idx] = sc.rv_histogram( (cdfY1E[1:] - cdfY1E[:-1],self.bins[idx]) )
	##}}}
	


