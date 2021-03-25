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
from SBCK.tools.__rv_extend import rv_mixture


###############
## Functions ##
###############

def gaussian_exp_mixture_1d( n_samples ):
	"""
	SBCK.datasets.gaussian_exp_mixtude_1d
	=====================================
	
	Build a test dataset such that:
	- X0 is normal distribution
	- X1 is a bimodal distribution, each mode is normal
	- Y0 is a mixture between an exponential and a normal distribution
	
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
	rvY0 = rv_mixture( [sc.expon(scale = 1),sc.norm(loc = 10,scale = 1)] , weights = [0.75,0.25] )
	rvX0 = sc.norm( loc = 7 , scale = 1 )
	rvX1 = rv_mixture( [sc.norm(loc = 5,scale = 1),sc.norm(loc = 9 , scale = 1)] , weights = [0.75,0.25] )
	Y0 = rvY0.rvs(n_samples).reshape(-1,1)
	X0 = rvX0.rvs(n_samples).reshape(-1,1)
	X1 = rvX1.rvs(n_samples).reshape(-1,1)
	
	return Y0,X0,X1


def gaussian_VS_exp_1d( n_samples ):
	"""
	SBCK.datasets.gaussian_VS_exp_1d
	================================
	
	Build a test dataset such that:
	- X0 is an exponential distribution
	- X1 is an exponential distribution
	- Y0 is a normal distribution
	
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
	
	Y0 = np.random.normal( loc = 5 , scale = 1 , size = (n_samples,1) )
	X0 = np.random.exponential( scale = 1   , size =  (n_samples,1) )
	X1 = np.random.exponential( scale = 0.5 , size =  (n_samples,1) )
	
	return Y0,X0,X1

