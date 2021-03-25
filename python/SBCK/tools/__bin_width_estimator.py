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


###############
## Functions ##
###############

def bin_width_estimator( X , method = "auto" ):
	"""
	SBCK.tools.bin_width_estimator
	==============================
	
	Estimate the width of the bin to build an histogram of X
	
	Parameters
	----------
	X    : np.array[ shape = (n_samples,n_features) ] or list(np.array)
		A dataset or a list of dataset X containing n_samples observations of n_features random variables.
	method : string = [ "auto" , "Sturges" , "FD" ]
		Method to estimate bin_width. If method == "auto", "Sturges" is selected if n_samples < 1000, else "FD"
	
	Returns
	-------
	bin_width : np.array[ shape = (n_features) ]
		bin_width of each features.
	"""
	
	if type(X) == list:
		return np.min( [ bin_width_estimator( x , method ) for x in X ] , axis = 0 )
	
	if X.ndim == 1 : X = X.reshape(-1,1)
	
	if method == "auto":
		method = "Sturges" if X.shape[0] < 1000 else "FD"
	
	if method == "Sturges":
		nh = np.log2( X.shape[0] ) + 1.
		bin_width = np.zeros(X.shape[1]) + 1. / nh
	else:
		bin_width = 2. * ( np.percentile( X , q = 75 , axis = 0 ) - np.percentile( X , q = 25 , axis = 0 ) ) / np.power( X.shape[0] , 1. / 3. )
	
	return bin_width

