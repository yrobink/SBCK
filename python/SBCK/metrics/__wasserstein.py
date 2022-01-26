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
from SBCK.tools.__OT import OTNetworkSimplex
from .__decorators import _to_SparseHist


##############
## Function ##
##############

@_to_SparseHist
def wasserstein( muX , muY , p = 2. , ot = OTNetworkSimplex() , metric = "euclidean" ):
	"""
	Description
	===========
	Compute the Wasserstein metric between two sparse histograms. If ot is a Sinkhorn algorithm, the dissimilarity is returned.
	
	Parameters
	----------
	muX      : SBCK.SparseHist or np.array
		Histogram or dataset
	muY      : SBCK.SparseHist or np.array
		Histogram or dataset
	p      : float
		Power of the cost function
	metric : str or callable
		See scipy.spatial.distance.pdist
	
	Return
	------
	cost   : float
		Wasserstein cost between muX and muY
	
	References
	----------
	[1] Wasserstein, L. N. (1969). Markov processes over denumerable products of spaces describing large systems of automata. Problems of Information Transmission, 5(3), 47-52.
	
	"""
	
	cost = lambda OT : np.sqrt(np.sum(OT.P * OT.C))
	
	ot.power = p 
	ot.fit( muX , muY )
	if type(ot) == OTNetworkSimplex:
		w = cost(ot)
		if not ot.state:
			w = np.nan
		if not abs(ot.P.sum() - 1) < 1e-6:
			w = np.nan
		return w
	else:
		costXY = cost(ot)
		ot.fit( muX , muX )
		costXX = cost(ot)
		ot.fit( muY , muY )
		costYY = cost(ot)
		return costXY - (costXX + costYY) / 2


