
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

class PPPFunctionLink(PrePostProcessing):##{{{
	
	def __init__( self , transform_ , itransform_ , *args , cols = None , **kwargs ):
		PrePostProcessing.__init__( self , *args , **kwargs )
		self._transform  = transform_
		self._itransform = itransform_
		self._cols = cols
		if cols is not None:
			self._cols = np.array( [cols] , dtype = int ).squeeze()
	
	def transform( self , X ):
		if self._cols is None:
			return self._transform(X)
		Xt = X
		Xt[:,self._cols] = self._transform(X[:,self._cols])
		return Xt
	
	def itransform( self , Xt ):
		if self._cols is None:
			return self._itransform(Xt)
		X = Xt
		X[:,self._cols] = self._itransform(Xt[:,self._cols])
		return Xt
##}}}

class PPPSquareLink(PPPFunctionLink):##{{{
	def __init__( self , *args , cols = None , **kwargs ):
		transform  = lambda x : x**2
		itransform = lambda x : np.where( x > 0 , np.sqrt(np.abs(x)) , - np.sqrt(np.abs(x)))
		PPPFunctionLink.__init__( self , transform , itransform , *args , cols = cols , **kwargs )
##}}}

class PPPLogLinLink(PPPFunctionLink):##{{{
	def __init__( self , *args , cols = None , **kwargs ):
		transform  = lambda x: np.where( (0 < x) & (x < 1) , np.log( np.where( x > 0 , x , np.nan ) ) , x - 1 )
		itransform = lambda x: np.where( x < 0 , np.exp(x) , x + 1 )
		PPPFunctionLink.__init__( self , transform , itransform , *args , cols = cols , **kwargs )
##}}}

