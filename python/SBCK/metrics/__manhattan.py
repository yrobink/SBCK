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

from .__minkowski import minkowski


##############
## Function ##
##############

def manhattan( muX , muY ):
	"""
	Description
	===========
	Manhattan distance between SparseHist
	
	Parameters
	----------
	muX      : SBCK.SparseHist or np.array
		Histogram or dataset
	muY      : SBCK.SparseHist or np.array
		Histogram or dataset

	Return
	------
	cost   : float
		Manhattan distance between muX and muY
	"""

	return minkowski( muX , muY , p = 1. )
