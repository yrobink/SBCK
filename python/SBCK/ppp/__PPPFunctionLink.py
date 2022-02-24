
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
	"""
	SBCK.ppp.PPPFunctionLink
	========================
	
	This class is used to define pre/post processing class with a link function
	and its inverse. See also the PrePostProcessing documentation
	
	>>> ## Start with data
	>>> Y0,X0,X1 = SBCK.datasets.like_tas_pr(2000)
	>>> 
	>>> ## Define the link function
	>>> transform  = lambda x : x**3
	>>> itransform = lamnda x : x**(1/3)
	>>> 
	>>> ## And the PPP method
	>>> ppp = SBCK.ppp.PPPFunctionLink( bc_method = SBCK.CDFt ,
	>>>                                transform_ = transform ,
	>>>                               itransform_ = itransform )
	>>> 
	>>> ## And now the correction
	>>> ## Bias correction
	>>> ppp.fit(Y0,X0,X1)
	>>> Z = ppp.predict(X1,X0)
	
	"""
	def __init__( self , transform_ , itransform_ , *args , cols = None , **kwargs ):
		"""
		Constructor
		===========
		
		Arguments
		---------
		transform_: [callable]
			Function to transform the data
		itransform_: [callable]
			Function to inverse the transform of the data
		cols: [int or array of int]
			The columns to apply the SSR
		isaved: str
			Choose the threshold used for inverse transform. Can be "Y0", "X0"
			or "X1"
		*args:
			All others arguments are passed to SBCK.ppp.PrePostProcessing
		*kwargs:
			All others arguments are passed to SBCK.ppp.PrePostProcessing
		"""
		PrePostProcessing.__init__( self , *args , **kwargs )
		self._transform  = transform_
		self._itransform = itransform_
		self._cols = cols
		if cols is not None:
			self._cols = np.array( [cols] , dtype = int ).squeeze()
	
	def transform( self , X ):
		"""
		Apply the transform
		"""
		if self._cols is None:
			return self._transform(X)
		Xt = X
		Xt[:,self._cols] = self._transform(X[:,self._cols])
		return Xt
	
	def itransform( self , Xt ):
		"""
		Apply the inverse transform
		"""
		if self._cols is None:
			return self._itransform(Xt)
		X = Xt
		X[:,self._cols] = self._itransform(Xt[:,self._cols])
		return Xt
##}}}

class PPPSquareLink(PPPFunctionLink):##{{{
	"""
	SBCK.ppp.PPPSquareLink
	======================
	
	Square link transform, i.e.:
	- transform is given by lambda x: x**2
	- inverse transform is given by lambda x: sign(x) * sqrt(abs(x))
	
	"""
	
	def __init__( self , *args , cols = None , **kwargs ):
		"""
		Constructor
		===========
		
		Arguments
		---------
		cols: [int or array of int]
			The columns to apply the Link function
		*args:
			All others arguments are passed to SBCK.ppp.PrePostProcessing
		*kwargs:
			All others arguments are passed to SBCK.ppp.PrePostProcessing
		"""
		transform  = lambda x : x**2
		itransform = lambda x : np.where( x > 0 , np.sqrt(np.abs(x)) , - np.sqrt(np.abs(x)))
		PPPFunctionLink.__init__( self , transform , itransform , *args , cols = cols , **kwargs )
##}}}

class PPPLogLinLink(PPPFunctionLink):##{{{
	"""
	SBCK.ppp.PPPLogLinLink
	======================
	
	Log linear link transform, i.e.:
	- transform is given by log(x) if 0 < x < 1, else x - 1
	- inverse transform is given by exp(x) if x < 0, else x + 1
	
	"""
	def __init__( self , *args , cols = None , **kwargs ):
		"""
		Constructor
		===========
		
		Arguments
		---------
		cols: [int or array of int]
			The columns to apply the Link function
		*args:
			All others arguments are passed to SBCK.ppp.PrePostProcessing
		*kwargs:
			All others arguments are passed to SBCK.ppp.PrePostProcessing
		"""
		transform  = lambda x: np.where( (0 < x) & (x < 1) , np.log( np.where( x > 0 , x , np.nan ) ) , x - 1 )
		itransform = lambda x: np.where( x < 0 , np.exp(x) , x + 1 )
		PPPFunctionLink.__init__( self , transform , itransform , *args , cols = cols , **kwargs )
##}}}

