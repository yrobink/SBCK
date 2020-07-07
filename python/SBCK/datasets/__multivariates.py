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

import numpy as np
import sklearn.datasets as skd

###############
## Functions ##
###############

def gaussian_exp_2d(n_samples):##{{{
	"""
	SBCK.datasets.gaussian_exp_2d
	=============================
	
	Build a bivariate test dataset.
	
	Parameters
	----------
	n_samples : integer
		Number of samples in X0, X1 and Y0
	
	Returns
	-------
	Y0,X0,X1 : tuple
		- Y0 reference dataset in calibration period, exp x norm
		- X0 biased dataset in calibration period, norm x exp
		- X1 biased dataset in projection period, (norm + 5) x exp
	"""
	X0 = np.hstack( ( np.random.normal( size = (n_samples,1) )           , np.random.exponential( size = (n_samples,1)  ) ) )
	Y0 = np.hstack( ( np.random.exponential( size = (n_samples,1)  )     , np.random.normal( size = (n_samples,1) ) ) )
	X1 = np.hstack( ( np.random.normal( size = (n_samples,1) , loc = 5 ) , np.random.exponential( size = (n_samples,1)  ) ) )
	
	return Y0,X0,X1
##}}}

def gaussian_L_2d( n_samples ):##{{{
	"""
	SBCK.datasets.gaussian_L_2d
	===========================
	
	Build a bivariate test dataset.
	
	Parameters
	----------
	n_samples : integer
		Number of samples in X0, X1 and Y0
	
	Returns
	-------
	Y0,X0,X1 : tuple
		- Y0 reference dataset in calibration period, form in "L"
		- X0 biased dataset in calibration period, gaussian
		- X1 biased dataset in projection period, gaussian
	"""
	## Construction of X0 (biased period 0), X1 (biased period 1) and Y0 (reference period 0)
	size0  = int(n_samples/2)
	size1  = n_samples - int(n_samples/4)
	
	## Just a gaussian for X0
	X0 = np.random.multivariate_normal( mean = [0.,0.] , cov = np.identity(2) , size = n_samples )
	
	## A lightly complex gaussian for X1
	X1 = np.random.multivariate_normal( mean = [1.,2.] , cov = [ [2.,0] , [0,0.5] ] , size = n_samples )
	
	## A very complex law for Y0
	Y0 = np.zeros( (n_samples,2) )
	Y0[:size0,:]      = np.random.multivariate_normal( mean = [7.,7.]   , cov = [[2,0],[0,0.5]]   , size = size0 )
	Y0[size0:size1,:] = np.random.multivariate_normal( mean = [5.,9.]   , cov = [[0.5,0],[0,2]]   , size = size1 - size0 )
	Y0[size1:]        = np.random.multivariate_normal( mean = [5.,12.5] , cov = [[0.2,0],[0,0.2]] , size = n_samples - size1 )
	meanY0 = np.mean( Y0 , axis = 0 )
	meanX0 = np.mean( X0 , axis = 0 )
	Y0 = np.apply_along_axis( lambda x : x - meanY0 + meanX0 , 1 , Y0 )
	
	return Y0,X0,X1
##}}}

def bimodal_reverse_2d( n_samples ):##{{{
	"""
	SBCK.datasets.bimodal_reverse_2d
	================================
	
	Build a test dataset such that:
	- X0 is a bimodal bivariate normal distribution, with different covariance matrix for each mode. Modes are close
	- X1 is a bimodal bivariate normal distribution, with different covariance matrix for each mode. Modes are differents
	- Y0 is a bimodal bivariate normal distribution, two modes are the same, but are orthogonal to X0 and X1
	
	Parameters
	----------
	n_samples : integer
		Number of samples in X0, X1 and Y0
	
	Returns
	-------
	Y0,X0,X1 : tuple
		- Y0 reference dataset in calibration period
		- X0 biased dataset in calibration period
		- X1 biased dataset in projection period
	"""
	drawn = int(n_samples/2)
	drawn = [drawn,n_samples-drawn]
	lmY0   = [ np.array([5,-3]) , np.array( [-3,3] ) ]
	lcovY0 = [ 0.9 * np.identity(2) , np.identity(2) ]
	lmX0   = [ np.zeros(2) , np.array( [2,2] ) ]
	lcovX0 = [ np.identity(2) , 0.5 * np.identity(2) ]
	lmX1   = [ np.zeros(2) - 1. , np.array( [5,5] ) ]
	lcovX1 = [ np.identity(2)  * 2 , 0.1 * np.identity(2) ]
	Y0     = np.vstack( [ np.random.multivariate_normal( mean = m , cov = cov , size = draw ) for m,cov,draw in zip(lmY0,lcovY0,drawn) ] )
	X0     = np.vstack( [ np.random.multivariate_normal( mean = m , cov = cov , size = draw ) for m,cov,draw in zip(lmX0,lcovX0,drawn) ] )
	X1     = np.vstack( [ np.random.multivariate_normal( mean = m , cov = cov , size = draw ) for m,cov,draw in zip(lmX1,lcovX1,drawn) ] )
	return Y0,X0,X1
##}}}

def gaussian_dd( n_samples , n_features = 2 ):##{{{
	"""
	SBCK.datasets.gaussian_dd
	=========================
	
	Build a test dataset such that X0, X1 and Y0 are multivariate normal distribution.
	
	Parameters
	----------
	n_samples : integer
		Number of samples in X0, X1 and Y0
	n_features : integer
		dimension, default is 2
	
	Returns
	-------
	Y0,X0,X1 : tuple
		- Y0 reference dataset in calibration period
		- X0 biased dataset in calibration period
		- X1 biased dataset in projection period
	"""
	X0 = np.random.multivariate_normal( mean = np.zeros(n_features)     , cov = skd.make_spd_matrix(n_features) , size = n_samples )
	X1 = np.random.multivariate_normal( mean = np.zeros(n_features) + 5 , cov = skd.make_spd_matrix(n_features) , size = n_samples )
	Y0 = np.random.multivariate_normal( mean = np.zeros(n_features) - 2 , cov = skd.make_spd_matrix(n_features) , size = n_samples )
	return Y0,X0,X1
##}}}

