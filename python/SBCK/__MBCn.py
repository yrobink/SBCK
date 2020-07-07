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
from .__QM import QM
from .__QDM import QDM
from .metrics.__wasserstein import wasserstein
from .tools.__SlopeStoppingCriteria import SlopeStoppingCriteria


###########
## Class ##
###########

class MBCn:
	"""
	SBCK.MBCn
	=========
	
	Description
	-----------
	MBCn Bias correction method, see [1]
	
	References
	----------
	[1] Cannon, Alex J.: Multivariate quantile mapping bias correction: an N-dimensional probability density function transform for climate model simulations of multiple variables, Climate Dynamics, nb. 1, vol. 50, p. 31-49, 10.1007/s00382-017-3580-6
	"""
	def __init__( self , bc = QDM , metric = wasserstein , stopping_criteria = SlopeStoppingCriteria , stopping_criteria_params = { "minit" : 20 , "maxit" : 100 , "tol" : 1e-3 } , **kwargs ): ##{{{
		"""
		Initialisation of MBCn.
		
		Parameters
		----------
		bc  : Bias correction method
			Non stationary bias correction method, default is QDM
		metric : callable
			Callable between two matrices, used as criteria to dermined when stopped iteration. Default is Wasserstein distance
		stopping_criteria: a class to determine when stop iteration
			See note
		stopping_criteria_params : dict
			Params of stopping_criteria
		kwargs : others named arguments
			Passed to bias correction method
		
		Note
		----
		MBCn method used an alternance of random rotation and quantile mapping to perform the multivariate bias correction.
		At each step, the metric is used in calibration period to determine if correction is close or not of correction.
		The class SlopeStoppingCriteria compute the slope of time series of callable. When the slope is lower than tol,
		or maxit is atteigned, iterations are stopped. At least minit is performed.
		Attributes
		----------
		"""
		self.iter_stop = stopping_criteria(**stopping_criteria_params)
		self.metric = metric
		self.bc = bc
		self.bc_params = kwargs
		self._lbc  = []
	##}}}
	
	@property
	def maxit(self):##{{{
		return self.iter_stop.maxit
	##}}}
	
	@property
	def nit(self):##{{{
		return self.iter_stop.nit
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
		
		if Y0.ndim == 1: Y0 = Y0.reshape(-1,1)
		if X0.ndim == 1: X0 = X0.reshape(-1,1)
		if X1.ndim == 1: X1 = X1.reshape(-1,1)
		
		n_features = Y0.shape[1]
		self.iter_stop.initialize()
		
		## Generate orthogonal matrices: SO(n_features)
		self.ortho_mat  = sc.special_ortho_group.rvs( n_features , self.maxit )
		
		## Tips for performance, inverse + ortho of next in one pass
		self.tips = np.zeros(self.ortho_mat.shape)
		for i in range(self.maxit-1):
			self.tips[i,:,:] = self.ortho_mat[i+1,:,:] @ np.linalg.inv(self.ortho_mat[i,:,:])
		self.tips[-1,:,:] = np.linalg.inv(self.ortho_mat[-1,:,:])
		
		## Loop
		Z0_o = np.transpose( self.ortho_mat[0,:,:] @ X0.T )
		Z1_o = np.transpose( self.ortho_mat[0,:,:] @ X1.T )
		
		for i in self.iter_stop:
			Y0_o = np.transpose( self.ortho_mat[i,:,:] @ Y0.T )
			
			bc = self.bc(**self.bc_params)
			bc.fit( Y0_o , Z0_o , Z1_o )
			Z1_o,Z0_o = bc.predict(Z1_o,Z0_o)
			
			self.iter_stop.append(self.metric(Z0_o,Y0_o))
			
			self._lbc.append(bc)
			
			Z0_o = np.transpose( self.tips[i,:,:] @ Z0_o.T )
			Z1_o = np.transpose( self.tips[i,:,:] @ Z1_o.T )
		
		Z0 = np.transpose( np.linalg.inv(self.ortho_mat[self.nit,:,:]) @ Z0_o.T )
		Z1 = np.transpose( np.linalg.inv(self.ortho_mat[self.nit,:,:]) @ Z1_o.T )
		
		self.ortho_mat = self.ortho_mat[:self.nit,:,:]
		self.tips = self.tips[:self.nit,:,:]
		self.tips[-1,:,:] = np.linalg.inv(self.ortho_mat[-1,:,:])
		
		bc = self.bc(**self.bc_params)
		bc.fit( Y0 , Z0 , Z1 )
		self._lbc.append(bc)
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
		if X0 is None:
			return self._predict_X1(X1)
		else:
			return self._predict_X1_X0(X1,X0)
		
	##}}}
	
	def _predict_X1( self , X1 ):##{{{
		if X1.ndim == 1: X1 = X1.reshape(-1,1)
		
		Z1_o = np.transpose( self.ortho_mat[0,:,:] @ X1.T ) if X1 is not None else None
		
		for i in range(self.nit):
			Z1_o = self._lbc[i].predict(Z1_o)
			Z1_o = np.transpose( self.tips[i,:,:] @ Z1_o.T )
		
		Z1 = self._lbc[-1].predict(Z1_o)
		return Z1
	##}}}
	
	def _predict_X1_X0( self , X1 , X0 ):##{{{
		if X0.ndim == 1: X0 = X0.reshape(-1,1)
		if X1.ndim == 1: X1 = X1.reshape(-1,1)
		
		Z0_o = np.transpose( self.ortho_mat[0,:,:] @ X0.T )
		Z1_o = np.transpose( self.ortho_mat[0,:,:] @ X1.T )
		
		for i in range(self.nit):
			Z1_o,Z0_o = self._lbc[i].predict(Z1_o,Z0_o)
			Z0_o = np.transpose( self.tips[i,:,:] @ Z0_o.T )
			Z1_o = np.transpose( self.tips[i,:,:] @ Z1_o.T )
		
		return self._lbc[-1].predict(Z1_o,Z0_o)
	##}}}


