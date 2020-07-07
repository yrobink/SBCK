# -*- coding: utf-8 -*-
##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the SBCK (Statistical    ##
## Bias Correction Kit). This library makes it possible to perform bias         ##
## correction with non parametric methods, and give some metrics between Sparse ##
## Histogram is high dimensions.                                                ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## SBCK (Statistical Bias Correction Kit). Cette librairie permet d'appliquer   ##
## une correction de biais avec des méthodes non paramétriques, et propose      ##
## diverses metrique entre Histograme Sparse en haute dimension.                ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################

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



