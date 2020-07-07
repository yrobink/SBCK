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
from .tools.__rv_extend import rv_histogram
from .__QM import QM


###########
## Class ##
###########


class MRec:
	"""
	SBCK.MRec
	=========
	
	Description
	-----------
	MRec Bias correction method, see [1]
	
	References
	----------
	[1] Bárdossy, A. and Pegram, G.: Multiscale spatial recorrelation of RCM precipitation to produce unbiased climate change scenarios over large areas and small, Water Resources Research, 48, 9502–, https://doi.org/10.1029/2011WR011524, 2012.
	"""
	
	def __init__( self , distY = None , distX = None ):##{{{
		"""
		Initialisation of MRec.
		
		Parameters
		----------
		distY: describe the distribution of reference, see QM.
			If None, rv_histogram is used
		distX: describe the distribution of biased dataset, see QM
			If None, rv_histogram is used
		
		"""
		self._qmX0 = None
		self._qmX1 = None
		self._qmY0 = None
		self._S_CY0g  = None
		self._Si_CX0g = None
		self._re_un_mat = None
		self.n_features = 0
		self._distY = distY
		self._distX = distX
	##}}}
	
	def fit( self , Y0 , X0 , X1 ):##{{{
		"""
		Fit the MBCn model
		
		Parameters
		----------
		Y0 : np.array[ shape = (n_samples,n_features) ]
			Reference dataset during period 0
		X0 : np.array[ shape = (n_samples,n_features) ]
			Biased dataset during period 0
		X1	: np.array[ shape = (n_samples,n_features) ]
			Biased dataset during period 1
		"""
		if Y0.ndim == 1 : Y0 = Y0.reshape(-1,1)
		if X0.ndim == 1 : X0 = X0.reshape(-1,1)
		if X1.ndim == 1 : X1 = X1.reshape(-1,1)
		self.n_features = Y0.shape[1]
		
		
		## Kind of variables
		if self._distY is None: self._distY = [rv_histogram for _ in range(self.n_features)]
		if self._distX is None: self._distX = [rv_histogram for _ in range(self.n_features)]
		
		## Transform into Gaussian data
		self._qmY0 = QM( distY0 = sc.norm(0,1) , distX0 = self._distY )
		self._qmX0 = QM( distY0 = sc.norm(0,1) , distX0 = self._distX )
		self._qmX1 = QM( distY0 = sc.norm(0,1) , distX0 = self._distX )
		self._qmY0.fit( None , Y0 )
		self._qmX0.fit( None , X0 )
		self._qmX1.fit( None , X1 )
		Y0g = self._qmY0.predict(Y0)
		X0g = self._qmX0.predict(X0)
		X1g = self._qmX1.predict(X1)
		
		## Correlation matrix
		CY0g = np.corrcoef( Y0g.T )
		CX0g = np.corrcoef( X0g.T )
		
		## Squareroot matrix
		a_CY0g,d_CY0g,_ = np.linalg.svd(CY0g)
		self._S_CY0g = a_CY0g @ np.diag(np.sqrt(d_CY0g)) @ a_CY0g.T
		
		a_CX0g,d_CX0g,_ = np.linalg.svd(CX0g)
		self._Si_CX0g = a_CX0g @ np.diag( np.power(d_CX0g,-0.5) ) @ a_CX0g.T
		
		## Decor-recor-relation
		self._re_un_mat = self._S_CY0g @ self._Si_CX0g
		X0_recor = np.transpose( self._re_un_mat @ X0g.T )
		X1_recor = np.transpose( self._re_un_mat @ X1g.T )
		
		## Final QM
		self._qmY0 = QM( distY0 = self._distY , distX0 = sc.norm )
		self._qmY0.fit( Y0 , X0_recor )
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
		"""
		if X0 is not None and X0.ndim == 1 : X0 = X0.reshape(-1,1)
		if X1.ndim == 1 : X1 = X1.reshape(-1,1)
		
		X1g = self._qmX1.predict(X1)
		X1_recor = np.transpose( self._re_un_mat @ X1g.T )
		Z1 = self._qmY0.predict(X1_recor)
		
		Z0 = None
		if X0 is not None:
			X0g = self._qmX0.predict(X0)
			X0_recor = np.transpose( self._re_un_mat @ X0g.T )
			Z0 = self._qmY0.predict(X0_recor)
			return Z1,Z0
		return Z1
	##}}}


