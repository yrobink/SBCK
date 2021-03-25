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

import numpy as np
from .__decorators import _to_SparseHist


##############
## Function ##
##############

@_to_SparseHist
def chebyshev( muX , muY ):
	"""
	Description
	===========
	Chebyshev distance between SparseHist, defines by
	dist = max_{ij} |x_i-y_i|

	Parameters
	----------
	muX      : SBCK.SparseHist or np.array
		Histogram or dataset
	muY      : SBCK.SparseHist or np.array
		Histogram or dataset

	Return
	------
	cost   : float
		Minkowski distance between muX and muY
	"""
	dist = 0
	indx = muY.argwhere( muX.c )
	indy = muX.argwhere( muY.c )

	## Common elements of muX in muY
	g = np.argwhere( indx > -1 ).ravel()
	for i in g:
		dist = max( dist , np.abs( muX.p[i] - muY.p[indx[i]] ) )
	
	## Elements of muX not in muY
	g = np.argwhere(indx == -1).ravel()
	for i in g:
		dist = max( dist , muX.p[i] )

	## Elements of muY not in muX
	g = np.argwhere(indy == -1).ravel()
	for i in g:
		dist = max( dist , muY.p[i] )

	return dist
