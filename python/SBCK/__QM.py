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


