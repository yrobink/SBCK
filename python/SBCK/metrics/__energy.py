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
import scipy.spatial.distance as ssd
from .__decorators import _to_SparseHist


##############
## Function ##
##############

@_to_SparseHist
def energy( muX , muY , p = 2. , metric = "euclidean" ):
	"""
	Description
	===========
	Energy distance between sparse histograms

	Parameters
	----------
	muX      : SBCK.SparseHist or np.array
		Histogram or dataset
	muY      : SBCK.SparseHist or np.array
		Histogram or dataset
	p        : float
		Power of the metric function
	metric   : str or callable
		See scipy.spatial.distance.pdist
	
	Return
	------
	distance : float
		Estimation of energy distance
	"""
	sizeX = muX.size
	sizeY = muY.size
	XY = np.power( ssd.cdist( muX.c , muY.c , metric = metric ) , p ) * np.dot( muX.p.reshape( (sizeX,1) ) , muY.p.reshape( (1,sizeY) ) )
	XX = ssd.squareform( np.power( ssd.pdist( muX.c , metric = metric ) , p ) ) * np.dot( muX.p.reshape( (sizeX,1) ) , muX.p.reshape( (1,sizeX) ) )
	YY = ssd.squareform( np.power( ssd.pdist( muY.c , metric = metric ) , p ) ) * np.dot( muY.p.reshape( (sizeY,1) ) , muY.p.reshape( (1,sizeY) ) )

	return np.power( 2 * np.sum(XY) - np.sum(XX) - np.sum(YY) , 1. / p )

