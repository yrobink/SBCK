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
##{{{

## Scientific libraries
##=====================

import numpy as np
import scipy.stats as sc
import SBCK as bc
import SBCK.datasets as bcd

## Plot libraries ##
##==================

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

#mpl.rcParams['font.size'] = 30
#plt.rc('text',usetex=True)
#plt.rcParams['text.latex.unicode'] = True

##}}}
###############

###############
## Fonctions ##
###############


def fig_univariate():##{{{
	
	## Construction of biased and reference dataset
	Y0,X0,X1 = bcd.gaussian_exp_mixture_1d(10000)
	
	## Bias correction with CDFt
	cdft = bc.CDFt()
	cdft.fit( Y0 , X0 , X1 )
	Z1,Z0 = cdft.predict( X1 , X0 )
	Z1[Z1<0] = 0
	
	## Bias correction with QDM
	qdm = bc.QDM( "multiplicative" )
	qdm.fit( Y0 , X0 , X1 )
	Z1qdm = qdm.predict(X1)
	
	## Random Variable
	bins = np.arange( -1 , 14 , 0.1 )
	rvY0 = sc.rv_histogram( np.histogram( Y0 , bins ) )
	rvX0 = sc.rv_histogram( np.histogram( X0 , bins ) )
	rvX1 = sc.rv_histogram( np.histogram( X1 , bins ) )
	rvZ0 = sc.rv_histogram( np.histogram( Z0 , bins ) )
	rvZ1 = sc.rv_histogram( np.histogram( Z1 , bins ) )
	rvZ1qdm = sc.rv_histogram( np.histogram( Z1qdm , bins ) )
	
	## Plot
	xmin = min( [T.min() for T in [Y0,X0,X1,Z0,Z1,Z1qdm]] )
	xmax = max( [T.max() for T in [Y0,X0,X1,Z0,Z1,Z1qdm]] )
	xd   = 0.05 * (xmax - xmin)
	ylim = (0,0.8)
	
	bins = np.linspace( xmin - xd , xmax + xd , 80 )
	fig_factor = 0.3
	fig = plt.figure( figsize = ( fig_factor * 30 , fig_factor * 20) )
	
	ax = fig.add_subplot( 2 , 3 , 1 )
	ax.hist( X0 , bins = bins , color = "red" , density = True , alpha = 0.5 )
	ax.set_ylim( ylim )
	ax.set_title( r"$X_0$" )
	
	ax = fig.add_subplot( 2 , 3 , 2 )
	ax.hist( Z0 , bins = bins , color = "green" , density = True , alpha = 0.5 , label = "QM" )
	ax.set_ylim( ylim )
	ax.set_title( r"$Z_0$" )
	ax.legend( loc = "upper right" )
	
	ax = fig.add_subplot( 2 , 3 , 3 )
	ax.hist( Y0 , bins = bins , color = "blue" , density = True , alpha = 0.5 )
	ax.set_ylim( ylim )
	ax.set_title( r"$Y_0$" )
	
	ax = fig.add_subplot( 2 , 3 , 4 )
	ax.hist( X1 , bins = bins , color = "red" , density = True , alpha = 0.5 )
	ax.set_ylim( ylim )
	ax.set_title( r"$X_1$" )
	
	ax = fig.add_subplot( 2 , 3 , 5 )
	ax.hist( Z1    , bins = bins , color = "green"  , density = True , alpha = 0.5 , label = "CDFt" )
	ax.hist( Z1qdm , bins = bins , color = "orange" , density = True , alpha = 0.5 , label = "QDM" )
	ax.set_ylim( ylim )
	ax.set_title( r"$Z_1$" )
	ax.legend( loc = "upper right" )
	
	ax = fig.add_subplot( 2 , 3 , 6 )
	ax.plot( bins , rvY0.cdf(bins)    , color = "blue"   , label = r"$Y_0$" )
	ax.plot( bins , rvX0.cdf(bins)    , color = "red"    , label = r"$X_0$" )
	ax.plot( bins , rvX1.cdf(bins)    , color = "red"    , label = r"$X_1$" , linestyle = "--" )
	ax.plot( bins , rvZ0.cdf(bins)    , color = "green"  , label = r"$Z_0$" )
	ax.plot( bins , rvZ1.cdf(bins)    , color = "green" , label  = r"$Z_1$ (CDFt)" , linestyle = "--" )
	ax.plot( bins , rvZ1qdm.cdf(bins) , color = "orange" , label = r"$Z_1$ (QDM)"  , linestyle = "--" )
	ax.legend( loc = "lower right" , fontsize = 8 , ncol = 1 )
	ax.set_title( "CDF" )
	
	fig.set_tight_layout(True)
	plt.savefig( "univariate.png" )
##}}}

def fig_multivariate():##{{{
	
	Y0,X0,X1 = bcd.gaussian_L_2d(10000)
	
	lbcm    = [ bc.dOTC                                          , bc.MBCn , bc.MRec ]
	lkwargs = [ {"bin_width" : [0.1,0.1] , "cov_factor" : "std"} , {} , {}]
	lname   = ["dOTC"                                            , "MBCn","MRec"]
	
	for bcm,kwargs,name in zip(lbcm,lkwargs,lname):
	
		## Bias correction
		met = bcm( **kwargs )
		met.fit( Y0 , X0 , X1 )
		Z1,Z0 = met.predict(X1,X0)
		
		xylim0 = -10
		xylim1 = 10
		X,Y = np.mgrid[xylim0:xylim1:100j, xylim0:xylim1:100j]
		pos = np.vstack( [X.ravel(),Y.ravel()] )
		cmapX = plt.cm.Reds
		cmapY = plt.cm.Blues
		cmapZ = plt.cm.Greens
		
		## Plot
		fig_factor = 0.3
		fig = plt.figure( figsize = ( fig_factor * 30 , fig_factor * 20) )
		
		ax = fig.add_subplot( 2 , 3 , 1 )
		kde = sc.gaussian_kde(X0.T)
		ax.imshow( np.rot90(kde(pos).reshape(X.shape)) , cmap = cmapX , extent = [xylim0,xylim1,xylim0,xylim1] )
		ax.set_title( r"$X_0$" )
		
		ax = fig.add_subplot( 2 , 3 , 2 )
		kde = sc.gaussian_kde(Z0.T)
		ax.imshow( np.rot90(kde(pos).reshape(X.shape)) , cmap = cmapZ , extent = [xylim0,xylim1,xylim0,xylim1] )
		ax.set_title( r"$Z_0$ (" + name + ")" )
		
		ax = fig.add_subplot( 2 , 3 , 3 )
		kde = sc.gaussian_kde(Y0.T)
		ax.imshow( np.rot90(kde(pos).reshape(X.shape)) , cmap = cmapY , extent = [xylim0,xylim1,xylim0,xylim1] )
		ax.set_title( r"$Y_0$" )
		
		ax = fig.add_subplot( 2 , 3 , 4 )
		kde = sc.gaussian_kde(X1.T)
		ax.imshow( np.rot90(kde(pos).reshape(X.shape)) , cmap = cmapX , extent = [xylim0,xylim1,xylim0,xylim1] )
		ax.set_title( r"$X_1$" )
		
		ax = fig.add_subplot( 2 , 3 , 5 )
		kde = sc.gaussian_kde(Z1.T)
		ax.imshow( np.rot90(kde(pos).reshape(X.shape)) , cmap = cmapZ , extent = [xylim0,xylim1,xylim0,xylim1] )
		ax.set_title( r"$Z_1$ (" + name + ")" )
		
		fig.set_tight_layout(True)
		plt.savefig( "multivariate_" + name + ".png" )
##}}}

def fig_sparsehist():##{{{
	X = np.random.multivariate_normal( mean = np.zeros(2) , cov = np.identity(2) , size = 1000000 )
	muX = bc.tools.SparseHist( X , bc.tools.bin_width_estimator(X) )
	
	fig = plt.figure( figsize = (8,8) )
	ax = fig.add_subplot( 1 , 1 , 1 , projection = "3d" )
	ax.scatter( muX.c[:,0] , muX.c[:,1] , muX.p , c = muX.p , cmap = plt.cm.inferno )
	ax.set_xlabel( r"$x$" )
	ax.set_ylabel( r"$y$" )
	ax.set_zlabel( r"$\mathrm{Density}$" )
	
	fig.set_tight_layout(True)
	plt.savefig("SparseHist.png")

##}}}


##########
## main ##
##########

if __name__ == "__main__":
	fig_univariate()
	fig_multivariate()
	
	print("Done")



