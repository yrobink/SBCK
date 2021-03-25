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
from .__chebyshev  import chebyshev
from .__decorators import _to_SparseHist


##############
## Function ##
##############

@_to_SparseHist
def minkowski( muX , muY , p ):
	"""
	Description
	===========
	Minkowski distance between SparseHist, defines by
	dist^p = sum_{ij} |x_i-y_i|^p

	Parameters
	----------
	muX      : SBCK.SparseHist or np.array
		Histogram or dataset
	muY      : SBCK.SparseHist or np.array
		Histogram or dataset
	p      : float or np.inf (for Chebyshev distance)
		Power of the distance. If p = 2, it is euclidean distance

	Return
	------
	cost   : float
		Minkowski distance between muX and muY
	"""
	if p == np.inf:
		return chebyshev(muX,muY)
	
	dist = 0
	indx = muY.argwhere( muX.c )
	indy = muX.argwhere( muY.c )

	## Common elements of muX in muY
	ii = np.argwhere( indx > -1 ).ravel()
	dist += np.sum( np.power( np.abs( muX.p[ii] - muY.p[indx[ii]] ) , p ) )
	
	## Elements of muX not in muY
	dist += np.sum( np.power( np.abs( muX.p[ np.argwhere(indx == -1).ravel() ] ) , p ) )
	
	## Elements of muY not in muX
	dist += np.sum( np.power( np.abs( muY.p[ np.argwhere(indy == -1).ravel() ] ) , p ) )
	
	return np.power( dist , 1. / p )
