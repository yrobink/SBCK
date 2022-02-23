
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

###########
## Class ##
###########

class PrePostProcessing:##{{{
	
	"""
	SBCK.ppp.PrePostProcessing
	==========================
	
	Base class to apply pre / post processing operations before / after a bias
	correction method. This class is equivalent to apply the identity function,
	and is used to be herited by others PrePostProcessing classes.
	
	"""
	
	
	def __init__( self , bc_method = None , bc_method_kwargs = {} , pipe = [] , pipe_kwargs = [] ):##{{{
		"""
		Constructor
		===========
		
		Arguments
		---------
		bc_method: [SBCK.<Bias Correction class]
			A bias correction method, optional if this class is given in 'pipe'
		bc_method_kwargs: [dict]
			Keyword arguments given to the constructor of bc_method
		pipe: [list of PrePostProcessing class]
			List of preprocessing class to apply to data before / after bias
			correction.
		pipe_kwargs: [list of dict]
			List of keyword arguments to pass to each PreProcessing class
		
		"""
		if not type(pipe) == list:
			pipe = [pipe]
		if not type(pipe_kwargs) == list:
			pipe_kwargs = [pipe_kwargs]
		
		self._pipe = [ p(**kwargs) for p,kwargs in zip(pipe,pipe_kwargs) ]
		if bc_method is not None:
			self._bc_method  = bc_method( **bc_method_kwargs )
		##}}}
	
	def transform( self , X ):##{{{
		"""
		Transformation to apply before the bias correction method
		"""
		return X
	##}}}
	
	def itransform( self , X ):##{{{
		"""
		Transformation to apply after the bias correction method
		"""
		return X
	##}}}
	
	def _pipe_transform( self , X ):##{{{
		if X is None:
			return None
		Xt = X.copy()
		
		for p in self._pipe[::-1]:
			Xt = p.transform(Xt)
		
		Xt = self.transform(Xt)
		
		return Xt
	##}}}
	
	def _pipe_itransform( self , Xt ):##{{{
		if Xt is None:
			return None
		X = Xt.copy()
		
		X = self.itransform(X)
		for p in self._pipe:
			X = p.itransform(X)
		
		return X
	##}}}
	
	def fit( self , Y0 , X0 , X1 = None ):##{{{
		"""
		Fit the bias correction method after the pre-processing.
		"""
		
		Y0t = self._pipe_transform(Y0)
		X0t = self._pipe_transform(X0)
		X1t = self._pipe_transform(X1)
		
		if X1 is None:
			self._bc_method.fit( Y0t , X0t )
		else:
			self._bc_method.fit( Y0t , X0t , X1t )
	##}}}
	
	def predict( self , X1 = None , X0 = None ):##{{{
		"""
		Predict the bias correction method after the pre-processing, then apply
		the post-processing operation.
		"""
		
		X0t = self._pipe_transform(X0)
		X1t = self._pipe_transform(X1)
		Z0t = None
		Z1t = None
		
		if X0 is None:
			Z1t = self._bc_method.predict(X1t)
			return self._pipe_itransform(Z1t)
		elif X1 is None:
			Z0t = self._bc_method.predict(X0t)
			return self._pipe_itransform(Z0t)
		else:
			Z1t,Z0t = self._bc_method.predict(X1t,Z0t)
			Z1 = self._pipe_itransform(Z1t)
			Z0 = self._pipe_itransform(Z0t)
			return Z1,Z0
	##}}}
	
##}}}

