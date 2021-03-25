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


##############
## Function ##
##############

def entropy( muX , muY = None ):
	"""
	Description
	===========
	Compute the entropy of muX. If muY is not None, compute the Kullback-Leibler divergence

	Parameters
	----------
	muX      : SBCK.SparseHist or np.array
		Histogram or dataset
	muY      : SBCK.SparseHist or np.array
		Histogram or dataset

	Return
	------
	real   : float
		Entropy or Kullback-Leilbler divergence
	"""
	
	p0 = muX.p
	p1 = np.ones_like(p0) if muY is None else muY.p
	
	return np.sum( p0 * np.log( p0 / p1 ) )

