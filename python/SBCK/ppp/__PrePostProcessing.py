
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
	
	This base class can be considered as the identity pre-post processing, and
	is used to be herited by others pre/post processing class. The key ideas are:
	- A PrePostProcessing based class contains a bias correction method, initalized
	  by the `bc_method` argument, always available for all herited class
	- The `pipe` keyword is a list of pre/post processing class, applied one after
	  the other.
	
	Try with an example, start with a dataset similar to tas/pr:
	>>> Y0,X0,X1 = SBCK.datasets.like_tas_pr(2000)
	
	The first column is Gaussian, but the second is an exponential law with a Dirac
	mass at 0, represented the 0 of precipitations. For a quantile mapping
	correction in the calibration period, we just apply:
	>>> qm = SBCK.QM()
	>>> qm.fit(Y0,X0)
	>>> Z0 = qm.predict(X0)
	
	Now, if we want to pre-post process with the SSR method (0 are replaced by
	random values between 0 (excluded) and the minimal non zero value), we write:
	>>> ppp = SBCK.ppp.PPPSSR( bc_method = SBCK.QM , cols = [2] )
	>>> ppp.fit(Y0,X0)
	>>> Z0 = ppp.predict(X0)
	
	The SSR approach is applied only on the second column (the precipitation), and
	the syntax is the same than for a simple bias correction method.
	
	Imagine now that we want to apply the SSR, and to ensure the positivity of CDFt
	for precipitation, we also want to use the LogLinLink pre-post processing
	method. This can be done with the following syntax:
	>>> ppp = SBCK.ppp.PPPLogLinLink( bc_method = SBCK.CDFt , cols = [2] ,
	>>>                               pipe = [SBCK.ppp.PPPSSR] ,
	>>>                               pipe_kwargs = [{"cols" : 2}] )
	>>> ppp.fit(Y0,X0,X1)
	>>> Z = ppp.predict(X1,X0)
	
	With this syntax, the pre processing operation is
	PPPLogLinLink.transform(PPPSSR.transform(data)) and post processing operation
	PPPSSR.itransform(PPPLogLinLink.itransform(bc_data)). So the formula can read
	from right to left (as the mathematical composition). Note it is equivalent
	to define:
	>>> ppp = SBCK.ppp.PrePostProcessing$new( bc_method = SBCK.CDFt,
	>>>                          pipe = [SBCK.ppp.PPPLogLinLink,SBCK.ppp.PPPSSR],
	>>>                          pipe_kwargs = [ {"cols":2} , {"cols":2} ] )
	
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
		
		self._kind = None
		
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
	
	def _pipe_transform( self , X , kind ):##{{{
		if X is None:
			return None
		Xt = X.copy()
		
		self._kind = kind
		for p in self._pipe[::-1]:
			Xt = p.transform(Xt)
		
		Xt = self.transform(Xt)
		
		return Xt
	##}}}
	
	def _pipe_itransform( self , Xt , kind ):##{{{
		if Xt is None:
			return None
		X = Xt.copy()
		
		self._kind = kind
		X = self.itransform(X)
		for p in self._pipe:
			X = p.itransform(X)
		
		return X
	##}}}
	
	def fit( self , Y0 , X0 , X1 = None ):##{{{
		"""
		Fit the bias correction method after the pre-processing.
		"""
		
		Y0t = self._pipe_transform( Y0 , "Y0" )
		X0t = self._pipe_transform( X0 , "X0" )
		X1t = self._pipe_transform( X1 , "X1" )
		
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
		
		X0t = self._pipe_transform( X0 , "X0" )
		X1t = self._pipe_transform( X1 , "X1" )
		Z0t = None
		Z1t = None
		
		if X0 is None:
			Z1t = self._bc_method.predict(X1t)
			return self._pipe_itransform( Z1t , "X1" )
		elif X1 is None:
			Z0t = self._bc_method.predict(X0t)
			return self._pipe_itransform( Z0t , "X0" )
		else:
			Z1t,Z0t = self._bc_method.predict(X1t,X0t)
			Z1 = self._pipe_itransform( Z1t , "X1" )
			Z0 = self._pipe_itransform( Z0t , "X0" )
			return Z1,Z0
	##}}}
	
##}}}

