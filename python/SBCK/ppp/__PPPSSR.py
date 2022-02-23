
## Copyright(c) 2022 Yoann Robin
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
from .__PrePostProcessing import PrePostProcessing


###########
## Class ##
###########

class PPPSSR(PrePostProcessing): ##{{{
	
	def __init__( self , *args , cols = None , isaved = "Y0" , **kwargs ): ##{{{
		PrePostProcessing.__init__( self , *args , **kwargs )
		self.Xn        = None
		self._cols     = cols
		self._isaved   = isaved
		self._icurrent = -1
		
		if cols is not None:
			self._cols = np.array( [cols] , dtype = int ).squeeze()
		
		if isaved == "Y0":
			self._isaved = 0
		elif isaved == "X0":
			self._isaved = 1
		elif isaved == "X1":
			self._isaved = 2
	##}}}
	
	def transform( self , X ):##{{{
		self._icurrent += 1
		
		if self._cols is None:
			self._cols = [i for i in range(X.shape[1])]
		cols = self._cols
		
		Xn = np.nanmin( np.where( X[:,cols] > 0 , X[:,cols] , np.nan ) , axis = 0 )
		if self._isaved == self._icurrent:
			self.Xn = Xn
		
		ncols = cols.size
		Xt = X.copy()
		Xt[:,cols] = np.where( (X[:,cols] > Xn).reshape(-1,ncols) , X[:,cols].reshape(-1,ncols) , np.random.uniform( low = Xn / 100 , high = Xn , size = (X.shape[0],cols.size) ) ).squeeze()
		
		return Xt
	##}}}
	
	def itransform( self , Xt ):##{{{
		
		X = Xt.copy()
		cols = self._cols
		X[:,cols] = np.where( Xt[:,cols] > self.Xn , Xt[:,cols] , 0 )
		
		return X
		##}}}
	
##}}}

