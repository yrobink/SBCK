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

def fig_qm_cdft():##{{{
	
	## Construction of biased and reference dataset
	size = 10000
	dsize = 7500
	
	X0 = np.random.normal( loc = 7 , scale = 1 , size = (size,1) )
	
	X1 = np.zeros( (size,1) )
	X1[:(size-dsize),:] = np.random.normal( loc = 5 , scale = 1 , size = (size-dsize,1) )
	X1[(size-dsize):,:] = np.random.normal( loc = 9 , scale = 1 , size = (dsize,1) )
	
	Y0 = np.zeros( (size,1) )
	Y0[:dsize,:] = np.random.exponential( scale = 1 , size =  (dsize,1) )
	Y0[dsize:,:] = np.random.normal( loc = 10 , scale = 1 , size = (size-dsize,1) )
	
	## Construction of corrector for period 0
	qm = bc.QM()
	qm.fit(Y0,X0)
	
	## Correction of period 0
	Z0 = qm.predict(X0)
	
	## Construction of corrector for period 1
	cdft = bc.CDFt()
	cdft.fit( Y0 , X0 , X1 )
	
	## Correction of period 1
	Z1 = cdft.predict(X1)
	
	## Random Variable
	bins = np.arange( -1 , 14 , 0.1 )
	rvY0 = sc.rv_histogram( np.histogram( Y0 , bins ) )
	rvX0 = sc.rv_histogram( np.histogram( X0 , bins ) )
	rvX1 = sc.rv_histogram( np.histogram( X1 , bins ) )
	rvUX0 = sc.rv_histogram( np.histogram( Z0 , bins ) )
	rvUX1 = sc.rv_histogram( np.histogram( Z1 , bins ) )
	
	## Plot
	bins = cdft.bins[0]
	fig_factor = 0.3
	fig = plt.figure( figsize = ( fig_factor * 30 , fig_factor * 20) )
	
	ax = fig.add_subplot( 2 , 3 , 1 )
	ax.hist( X0 , bins = bins , color = "red" , density = True , alpha = 0.5 )
	ax.set_ylim( (0,0.8) )
	ax.set_title( r"$X_0$" )
	
	ax = fig.add_subplot( 2 , 3 , 2 )
	ax.hist( Z0 , bins = bins , color = "green" , density = True , alpha = 0.5 )
	ax.set_ylim( (0,0.8) )
	ax.set_title( r"$Z_0$" )
	
	ax = fig.add_subplot( 2 , 3 , 3 )
	ax.hist( Y0 , bins = bins , color = "blue" , density = True , alpha = 0.5 )
	ax.set_ylim( (0,0.8) )
	ax.set_title( r"$Y_0$" )
	
	ax = fig.add_subplot( 2 , 3 , 4 )
	ax.hist( X1 , bins = bins , color = "red" , density = True , alpha = 0.5 )
	ax.set_ylim( (0,0.8) )
	ax.set_title( r"$X_1$" )
	
	ax = fig.add_subplot( 2 , 3 , 5 )
	ax.hist( Z1 , bins = bins , color = "green" , density = True , alpha = 0.5 )
	ax.set_ylim( (0,0.8) )
	ax.set_title( r"$Z_1$" )
	
	ax = fig.add_subplot( 2 , 3 , 6 )
	ax.plot( bins , rvY0.cdf(bins)  , color = "blue"  , label = r"$Y_0$" )
	ax.plot( bins , rvX0.cdf(bins)  , color = "red"   , label = r"$X_0$" )
	ax.plot( bins , rvX1.cdf(bins)  , color = "red"   , linestyle = "--" , label = r"$X_1$" )
	ax.plot( bins , rvUX0.cdf(bins) , color = "green" , label = r"$Z_0$" )
	ax.plot( bins , rvUX1.cdf(bins) , color = "green" , linestyle = "--" , label = r"$Z_1$" )
	ax.legend( loc = "upper left" )
	ax.set_title( "CDF" )
	
	fig.set_tight_layout(True)
	plt.savefig( "qm_cdft.png" )
##}}}

def fig_otc_dotc():##{{{
	## Construction of X0 (biased period 0), X1 (biased period 1) and Y0 (reference period 0)
	factor = 10
	size   = factor * 2000
	sized0 = factor * 1000
	sized1 = factor * 1500
	
	## Just a gaussian for X0
	X0 = np.random.multivariate_normal( mean = [0.,0.] , cov = np.identity(2) , size = size )
	Y = np.zeros( (size,2) )
	
	## A lightly complex gaussian for X1
	X1 = np.random.multivariate_normal( mean = [1.,2.] , cov = [ [2.,0] , [0,0.5] ] , size = size )
	
	## A very complex law for Y0
	Y[:sized0,:] = np.random.multivariate_normal( mean = [7.,7.] , cov = np.array( [2,0,0,0.5] ).reshape( (2,2) ) , size = sized0 )
	Y[sized0:sized1,:] = np.random.multivariate_normal( mean = [5.,9.] , cov = np.array( [0.5,0,0,2] ).reshape( (2,2) ) , size = (sized1-sized0) )
	Y[sized1:] = np.random.multivariate_normal( mean = [5.,12.5] , cov = 0.2 * np.identity(2) , size = (size-sized1) )
	meanY = np.mean( Y , axis = 0 )
	meanX = np.mean( X0 , axis = 0 )
	Y = np.apply_along_axis( lambda x : x - meanY + meanX , 1 , Y )
	
	## Construction of corrector period0
	otc = bc.OTC( bin_width = [ 0.1 , 0.1 ] )
	otc.fit( Y , X0 )
	
	## Correction period0
	Z0 = otc.predict(X0)
	
	## Construction of corrector period1
	cov_factor = "std"
#	cov_factor = "cholesky"
#	cov_factor = "identity"
	dotc = bc.dOTC( bin_width = [0.1,0.1] , cov_factor = cov_factor )
	dotc.fit( Y , X0 , X1 )
	
	## Correction period1
	Z1 = dotc.predict(X1)
	
	## Pearson correlation
	pY,_   = sc.spearmanr( Y[:,0] , Y[:,1] )
	pX0,_  = sc.spearmanr( X0[:,0] , X0[:,1] )
	pX1,_  = sc.spearmanr( X1[:,0] , X1[:,1] )
	pUX0,_ = sc.spearmanr( Z0[:,0] , Z0[:,1] )
	pUX1,_ = sc.spearmanr( Z1[:,0] , Z1[:,1] )
	
	## Histogram for plot
	bins = [ np.arange( -6 , 10 , 0.1 ) for i in range(2) ]
	extent = [-8,8,-8,8]
	
	HX0,_,_ = np.histogram2d( X0[:,0] , X0[:,1] , bins = bins )
	HX0 = HX0 / np.sum(HX0)
	HX0[HX0 == 0] = np.nan
	
	HX1,_,_ = np.histogram2d( X1[:,0] , X1[:,1] , bins = bins )
	HX1 = HX1 / np.sum(HX1)
	HX1[HX1 == 0] = np.nan
	
	HY,_,_ = np.histogram2d( Y[:,0] , Y[:,1] , bins = bins )
	HY = HY / np.sum(HY)
	HY[HY == 0] = np.nan
	
	HZ0,_,_ = np.histogram2d( Z0[:,0] , Z0[:,1] , bins = bins )
	HZ0 = HZ0 / np.sum(HZ0)
	HZ0[HZ0 == 0] = np.nan
	
	HZ1,_,_ = np.histogram2d( Z1[:,0] , Z1[:,1] , bins = bins )
	HZ1 = HZ1 / np.sum(HZ1)
	HZ1[HZ1 == 0] = np.nan
	
	vmin = min( [ np.nanmin(X) for X in [ HX0 , HY , HZ0 ] ] )
	vmax = max( [ np.nanmax(X) for X in [ HX0 , HY , HZ0 ] ] )
	
	## Plot
	fig_factor = 0.3
	fig = plt.figure( figsize = ( fig_factor * 30 , fig_factor * 20) )
	
	ax = fig.add_subplot( 2 , 3 , 1 )
	ax.imshow( np.rot90(HX0) , cmap = plt.cm.inferno , extent = extent , vmin = vmin , vmax = vmax )
	ax.set_title( r"$X_0$" )
	
	ax = fig.add_subplot( 2 , 3 , 2 )
	ax.imshow( np.rot90(HZ0) , cmap = plt.cm.inferno , extent = extent , vmin = vmin , vmax = vmax )
	ax.set_title( r"$Z_0$" )
	
	ax = fig.add_subplot( 2 , 3 , 3 )
	ax.imshow( np.rot90(HY) , cmap = plt.cm.inferno , extent = extent , vmin = vmin , vmax = vmax )
	ax.set_title( r"$Y_0$" )
	
	ax = fig.add_subplot( 2 , 3 , 4 )
	ax.imshow( np.rot90(HX1) , cmap = plt.cm.inferno , extent = extent , vmin = vmin , vmax = vmax )
	ax.set_title( r"$X_1$" )
	
	ax = fig.add_subplot( 2 , 3 , 5 )
	ax.imshow( np.rot90(HZ1) , cmap = plt.cm.inferno , extent = extent , vmin = vmin , vmax = vmax )
	ax.set_title( r"$Z_1$" )
	
	fig.set_tight_layout(True)
	plt.savefig( "otc_dotc.png" )
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
	fig_qm_cdft()
	fig_otc_dotc()
	fig_sparsehist()
	
	print("Done")



