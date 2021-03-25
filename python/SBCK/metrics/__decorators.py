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
from SBCK.tools.__bin_width_estimator import bin_width_estimator
from SBCK.tools.__tools_cpp import SparseHist



##############
## Function ##
##############

def _to_SparseHist( func ):
	def wrapper( *args , **kwargs ):
		muX = args[0]
		muY = args[1]
		testX = type(muX) == SparseHist
		testY = type(muY) == SparseHist
		
		if not testX and testY:
			muXX = SparseHist( muX , muY.bin_width )
			return func( muXX , muY , **kwargs )
		elif testX and not testY:
			muYY = SparseHist( muY , muX.bin_width )
			return func( muX , muYY , **kwargs )
		elif not testX and not testY:
			bw = np.array( [bin_width_estimator([muX,muY])] )
			bw = bw.reshape(bw.size)
			muXX = SparseHist( muX , bw )
			muYY = SparseHist( muY , bw )
			return func( muXX , muYY , **kwargs )
		
		return func( muX , muY , **kwargs )
	return wrapper


