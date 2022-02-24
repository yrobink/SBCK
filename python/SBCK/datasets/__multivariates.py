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

import numpy as np
import scipy.stats as sc
import sklearn.datasets as skd

from ..__QM import QM


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

def like_tas_pr( n_samples ):##{{{
	"""
	SBCK.datasets.like_tas_pr
	=========================
	
	Build a test dataset such that X0, X1 and Y0 are similar to temperature/
	precipitation.
	
	The method is the following:
	- Data from a multivariate normal law (dim = 2) are drawn
	- The quantile mapping is used to map the last column into the exponential law
	- Values lower than a fixed quantile are replaced by 0
	
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
	n_dim = 2
	
	mX0   = np.array([5 for _ in range(n_dim-1)] + [0])
	mX1   = np.array([8 for _ in range(n_dim-1)] + [0])
	mY0   = np.zeros(n_dim)
	covX0 = skd.make_spd_matrix(n_dim)
	covX1 = skd.make_spd_matrix(n_dim)
	covY0 = skd.make_spd_matrix(n_dim)
	
	X0 = np.random.multivariate_normal( mean = mX0 , cov = covX0 , size = n_samples )
	X1 = np.random.multivariate_normal( mean = mX1 , cov = covX1 , size = n_samples )
	Y0 = np.random.multivariate_normal( mean = mY0 , cov = covY0 , size = n_samples )
	
	qm = QM( distY0 = sc.expon( scale = 1 ) )
	qm.fit(None,Y0[:,-1])
	Y0[:,-1] = qm.predict(Y0[:,-1]).squeeze()
	
	qm = QM( distY0 = sc.expon( scale = 0.5 ) )
	qm.fit(None,X0[:,-1])
	X0[:,-1] = qm.predict(X0[:,-1]).squeeze()
	
	qm = QM( distY0 = sc.expon( scale = 1 ) )
	qm.fit(None,X1[:,-1])
	X1[:,-1] = qm.predict(X1[:,-1]).squeeze()
	
	X0[X0[:,-1] < np.quantile(X0[:,-1],0.05),-1] = 0
	X1[X0[:,-1] < np.quantile(X1[:,-1],0.10),-1] = 0
	Y0[Y0[:,-1] < np.quantile(Y0[:,-1],0.35),-1] = 0
	
	return Y0,X0,X1
##}}}

