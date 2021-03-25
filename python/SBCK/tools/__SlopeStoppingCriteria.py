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

import numpy       as np
import scipy.stats as sc

###########
## Class ##
###########

class SlopeStoppingCriteria:
	def __init__( self , minit , maxit , tol ):
		self.minit    = minit
		self.maxit    = maxit
		self.nit      = -1
		self.tol      = tol
		self.stop     = False
		self.criteria = list()
		self.slope    = list()
	
	def initialize(self):
		self.nit      = -1
		self.stop     = False
		self.criteria = list()
		self.slope    = list()
	
	def append( self , value ):
		self.criteria.append(value)
		if self.nit > self.minit:
			slope,_,_,_,_ = sc.linregress( range(len(self.criteria)) , self.criteria )
			self.stop = np.abs(slope) < self.tol
			self.slope.append(slope)
	
	def __iter__(self):
		return self
	
	def __next__(self):
		self.nit += 1
		if not self.nit < self.maxit-1:
			self.stop = True
		if not self.stop:
			return self.nit
		raise StopIteration
