
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
###############

import sys
import numpy as np
import scipy.stats as sc
import statsmodels.tsa.stattools as stt

import SBCK as bc
import SBCK.tools as bct
import SBCK.metrics as bcm
import SBCK.datasets as bcd

import matplotlib as mpl
import matplotlib.pyplot as plt


#############
## Classes ##
#############

class AR:##{{{
	"""
	Class generating an Auto Regressive process (AR process), i.e.:
	x(n) = loc + phi[0]x(n-1) + phi[1] x(n-2) + ... + phi[k-1] x(n-k)
	
	"""
	def __init__( self , phi , loc = 0 , sigma = 1 , init = None ):##{{{
		"""
		Constructor of AR class
		
		phi : np.array
			Coefficients of AR process
		loc : float
			Center parameters
		sigma : float
			scale of white noise
		"""
		self.loc = loc
		self.phi = np.array(phi)
		self.noise = lambda : np.random.normal( size = 1 , scale = sigma )
		self.memory = self._init_memory() if init is None else init
		self._optim = None
	##}}}
	
	def _init_memory( self ):##{{{
		memory = np.zeros_like(self.phi)
		memory[0] = self.noise()
		for i in range(1,self.phi.size):
			memory[i] = np.sum( memory[:i] * self.phi[:i] ) + self.noise()
		return memory
	##}}}
	
	## Methods overloaded
	def __call__( self ):##{{{
		"""
		Return next value of process
		"""
		Xn = self.loc + np.sum( self.memory * self.phi ) + self.noise()
		self.memory[1:] = self.memory[:-1]
		self.memory[0] = Xn
		return Xn
	##}}}
	
	def plot_pacf( lX , color , ax = None , nlags = 7 , alpha = 0.1 ):##{{{
		
		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(1,1,1)
		
		##
		n_X = len(lX)
		lpac = []
		for X in lX:
			lpac.append( stt.pacf( X , nlags = nlags , alpha = alpha )[0] )
		
		lpac = np.array(lpac)[:,1:]
		
		for i in range(nlags):
			
			x_plt = np.linspace( i , i + 1 , n_X + 2 )[1:-1]
			for j in range(n_X):
				ax.vlines( x_plt[j] , 0 , lpac[j,i] , color = color[j] )
			ax.vlines( i , -1 , 1 , color = "black" , alpha = 0.5 )
		
		ax.fill_between( [0,nlags+1] , - alpha / 2 , alpha / 2 , color = "blue" , alpha = 0.3 )
		
		ax.hlines( 0 , 0 , nlags , color = "black" , alpha = 0.5 )
		ax.set_xlim( (0,nlags) )
		ax.set_ylim( (-1,1) )
		ax.set_xticks( np.arange(0,nlags,1) + 0.5 )
		ax.set_xticklabels( np.arange(1,nlags+1,1,dtype=np.int) )
	##}}}
##}}}


###############
## Functions ##
###############

def test_qm( plot = True ):##{{{
	## Dataset
	Y0,X0,X1 = bcd.gaussian_exp_mixture_1d(2000)
	
	## QM
	qm = bc.QM()
	qm.fit(Y0,X0)
	Z0 = qm.predict(X0)
	
	## CDFt
	cdft = bc.CDFt()
	cdft.fit( Y0 , X0 , X1 )
	Z1 = cdft.predict(X1)
	Z  = cdft.predict( X1 , X0 )
	
	
	## Plot
	if plot:
		
		xmin = min( [T.min() for T in [X0,X1,Y0,Z0,Z1]] )
		xmax = max( [T.max() for T in [X0,X1,Y0,Z0,Z1]] )
		bins = np.linspace( xmin , xmax , 80 )
		
		nrow,ncol,fs = 1,2,7
		fig = plt.figure( figsize = (fs*ncol,0.7*fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.hist( X0 , bins = bins , color = "red"   , alpha = 0.5 , density = True )
		ax.hist( Y0 , bins = bins , color = "blue"  , alpha = 0.5 , density = True )
		ax.hist( Z0 , bins = bins , color = "green" , alpha = 0.5 , density = True )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		ax.hist( X1 , bins = bins , color = "red"   , alpha = 0.5 , density = True )
		ax.hist( Z1 , bins = bins , color = "green" , alpha = 0.5 , density = True )
		
		plt.tight_layout()
		plt.show()
##}}}

def test_otc_univ( plot = True ):##{{{
	## Dataset
	Y0,X0,X1 = bcd.gaussian_exp_mixture_1d(2000)
	
	## BC
	otc = bc.OTC()
	otc.fit( Y0 , X0 )
	Z0 = otc.predict( X0 )
	
	dotc = bc.dOTC()
	dotc.fit( Y0 , X0 , X1 )
	Z1 = dotc.predict( X1 )
	Z = dotc.predict( X1 , X0 )
	
	
	## Plot
	if plot:
		xmin = min( [T.min() for T in [X0,X1,Y0,Z0,Z1]] )
		xmax = max( [T.max() for T in [X0,X1,Y0,Z0,Z1]] )
		bins = np.linspace( xmin , xmax , 50 )
		
		nrow,ncol,fs = 1,2,7
		fig = plt.figure( figsize = (fs*ncol,0.7*fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.hist( X0 , bins = bins , color = "red"   , alpha = 0.5 , density = True )
		ax.hist( Y0 , bins = bins , color = "blue"  , alpha = 0.5 , density = True )
		ax.hist( Z0 , bins = bins , color = "green" , alpha = 0.5 , density = True )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		ax.hist( X1 , bins = bins , color = "red"   , alpha = 0.5 , density = True )
		ax.hist( Z1 , bins = bins , color = "green" , alpha = 0.5 , density = True )
		
		plt.tight_layout()
		plt.show()
##}}}

def test_otc_biv( plot = True ):##{{{
	## Data
	Y0,X0,X1 = bcd.bimodal_reverse_2d(4000)
	
	## BC
	otc = bc.OTC()
	otc.fit( Y0 , X0 )
	Z0 = otc.predict( X0 )
	
	dotc = bc.dOTC()
	dotc.fit( Y0 , X0 , X1 )
	Z1 = dotc.predict( X1 )
	Z = dotc.predict( X1 , X0 )
	
	
	if plot:
		## Histogram
		xymin = np.min( [ Z.min() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		xymax = np.max( [ Z.max() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		bins = [ np.linspace( xymin , xymax , 100 ) for _ in range(2) ]
		
		lH = list()
		for K in [Y0,X0,X1,Z0,Z1]:
			H,_,_ = np.histogram2d( K[:,0] , K[:,1] , bins )
			H /= np.sum(H)
			H[H == 0] = np.nan
			vmax = np.nanmax(H)
			lH.append(H)
		
		
		## Plot
		cmapR = plt.cm.Reds
		cmapG = plt.cm.inferno
		cmapB = plt.cm.Blues
		nrow,ncol,fs = 2,2,5
		fig = plt.figure( figsize = (fs*ncol,fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.imshow( np.rot90(lH[0]) , cmap = cmapB , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		ax.imshow( np.rot90(lH[1]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		ax.imshow( np.rot90(lH[2]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 3 )
		ax.imshow( np.rot90(lH[3]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 4 )
		ax.imshow( np.rot90(lH[4]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		plt.show()
##}}}

def test_ECBC( plot = True ):##{{{
	
	## Data
	Y0,X0,X1 = bcd.bimodal_reverse_2d(4000)
	
	
	## Bias Correction
	irefs = [0]
	
	ecbc = bc.ECBC()
	ecbc.fit( Y0 , X0 , X1 )
	Z1,Z0 = ecbc.predict(X1,X0)
	
	
	## Histogram
	if plot:
		xymin = np.min( [ Z.min() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		xymax = np.max( [ Z.max() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		bins = [ np.linspace( xymin , xymax , 100 ) for _ in range(2) ]
		
		lH = list()
		for K in [Y0,X0,X1,Z0,Z1]:
			H,_,_ = np.histogram2d( K[:,0] , K[:,1] , bins )
			H /= np.sum(H)
			H[H == 0] = np.nan
			vmax = np.nanmax(H)
			lH.append(H)
		
		
		## Plot
		cmapR = plt.cm.Reds
		cmapG = plt.cm.inferno
		cmapB = plt.cm.Blues
		nrow,ncol,fs = 2,2,5
		fig = plt.figure( figsize = (fs*ncol,fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.imshow( np.rot90(lH[0]) , cmap = cmapB , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		ax.imshow( np.rot90(lH[1]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		ax.imshow( np.rot90(lH[2]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 3 )
		ax.imshow( np.rot90(lH[3]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 4 )
		ax.imshow( np.rot90(lH[4]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		plt.show()
##}}}

def test_qmrs( plot = True ):##{{{
	
	## Data
	Y0,X0,X1 = bcd.bimodal_reverse_2d(4000)
	
	
	## Bias Correction
	irefs = [0]
	
	qmrs = bc.QMrs( irefs = irefs )
	qmrs.fit( Y0 , X0 )
	Z0 = qmrs.predict(X0)
	
	r2d2 = bc.R2D2( irefs = irefs )
	r2d2.fit( Y0 , X0 , X1 )
	Z1 = r2d2.predict(X1)
	Z  = r2d2.predict(X1,X0)
	
	
	## Histogram
	if plot:
		xymin = np.min( [ Z.min() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		xymax = np.max( [ Z.max() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		bins = [ np.linspace( xymin , xymax , 100 ) for _ in range(2) ]
		
		lH = list()
		for K in [Y0,X0,X1,Z0,Z1]:
			H,_,_ = np.histogram2d( K[:,0] , K[:,1] , bins )
			H /= np.sum(H)
			H[H == 0] = np.nan
			vmax = np.nanmax(H)
			lH.append(H)
		
		
		## Plot
		cmapR = plt.cm.Reds
		cmapG = plt.cm.inferno
		cmapB = plt.cm.Blues
		nrow,ncol,fs = 2,2,5
		fig = plt.figure( figsize = (fs*ncol,fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.imshow( np.rot90(lH[0]) , cmap = cmapB , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		ax.imshow( np.rot90(lH[1]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		ax.imshow( np.rot90(lH[2]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 3 )
		ax.imshow( np.rot90(lH[3]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 4 )
		ax.imshow( np.rot90(lH[4]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		plt.show()
##}}}

def test_QDM( plot = True ):##{{{
	## Data
	Y0,X0,X1 = bcd.bimodal_reverse_2d(4000)
	
	## Bias Correction
	qdm = bc.QDM()
	qdm.fit( Y0 , X0 , X1 )
	Z1,Z0 = qdm.predict(X1,X0)
	
	
	if plot:
		## Histogram
		xymin = np.min( [ Z.min() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		xymax = np.max( [ Z.max() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		bins = [ np.linspace( xymin , xymax , 100 ) for _ in range(2) ]
		
		lH = list()
		for K in [Y0,X0,X1,Z0,Z1]:
			H,_,_ = np.histogram2d( K[:,0] , K[:,1] , bins )
			H /= np.sum(H)
			H[H == 0] = np.nan
			vmax = np.nanmax(H)
			lH.append(H)
		
		
		## Plot
		cmapR = plt.cm.Reds
		cmapG = plt.cm.inferno
		cmapB = plt.cm.Blues
		nrow,ncol,fs = 2,2,5
		fig = plt.figure( figsize = (fs*ncol,fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.imshow( np.rot90(lH[0]) , cmap = cmapB , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		ax.imshow( np.rot90(lH[1]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		ax.imshow( np.rot90(lH[2]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 3 )
		ax.imshow( np.rot90(lH[3]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 4 )
		ax.imshow( np.rot90(lH[4]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		plt.show()
##}}}

def test_MBCn( plot = True ):##{{{
	## Data
	Y0,X0,X1 = bcd.bimodal_reverse_2d(4000)
	
	## Bias Correction
	mbcn = bc.MBCn()
	mbcn.fit( Y0 , X0 , X1 )
	Z1,Z0 = mbcn.predict(X1,X0)
	
	if plot:
		## Histogram
		xymin = np.min( [ Z.min() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		xymax = np.max( [ Z.max() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		bins = [ np.linspace( xymin , xymax , 100 ) for _ in range(2) ]
		
		lH = list()
		for K in [Y0,X0,X1,Z0,Z1]:
			H,_,_ = np.histogram2d( K[:,0] , K[:,1] , bins )
			H /= np.sum(H)
			H[H == 0] = np.nan
			vmax = np.nanmax(H)
			lH.append(H)
		
		
		## Plot
		cmapR = plt.cm.Reds
		cmapG = plt.cm.inferno
		cmapB = plt.cm.Blues
		nrow,ncol,fs = 2,2,5
		fig = plt.figure( figsize = (fs*ncol,fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.imshow( np.rot90(lH[0]) , cmap = cmapB , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		ax.imshow( np.rot90(lH[1]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		ax.imshow( np.rot90(lH[2]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 3 )
		ax.imshow( np.rot90(lH[3]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 4 )
		ax.imshow( np.rot90(lH[4]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		plt.show()
##}}}

def test_MRec( plot = True ):##{{{
	## Data
	Y0,X0,X1 = bcd.bimodal_reverse_2d(4000)
	
	## Bias Correction
	mbcn = bc.MRec()
	mbcn.fit( Y0 , X0 , X1 )
	Z1 = mbcn.predict(X1)
	_,Z0 = mbcn.predict(X1,X0)
	
	if plot:
		## Histogram
		xymin = np.min( [ Z.min() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		xymax = np.max( [ Z.max() for Z in [Y0,X0,X1,Z0,Z1] ] ) 
		bins = [ np.linspace( xymin , xymax , 100 ) for _ in range(2) ]
		
		lH = list()
		for K in [Y0,X0,X1,Z0,Z1]:
			H,_,_ = np.histogram2d( K[:,0] , K[:,1] , bins )
			H /= np.sum(H)
			H[H == 0] = np.nan
			vmax = np.nanmax(H)
			lH.append(H)
		
		
		## Plot
		cmapR = plt.cm.Reds
		cmapG = plt.cm.inferno
		cmapB = plt.cm.Blues
		nrow,ncol,fs = 2,2,5
		fig = plt.figure( figsize = (fs*ncol,fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.imshow( np.rot90(lH[0]) , cmap = cmapB , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		ax.imshow( np.rot90(lH[1]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		ax.imshow( np.rot90(lH[2]) , cmap = cmapR , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 3 )
		ax.imshow( np.rot90(lH[3]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		ax = fig.add_subplot( nrow , ncol , 4 )
		ax.imshow( np.rot90(lH[4]) , cmap = cmapG , vmin = 0  , alpha = 1 , extent = [xymin,xymax,xymin,xymax] )
		
		plt.show()
##}}}

def test_tsbc( plot = True ):##{{{
	size = 2000
	
	## Generate AR3 processes
	##=======================
	ar3X = AR( loc = 0.2 , phi = [0.6,-0.2,0.1] , sigma = 1 )
	X  = np.array( [ ar3X() for _ in range(size) ] )
	
	ar3Y = AR( loc = 0 , phi = [-0.3,0.4,-0.2] , sigma = 0.7 )
	Y  = np.array( [ ar3Y() for _ in range(size) ] ) + 5
	
	
	## Bias correction
	##================
	nlags = 50
	
	tsbc = bc.TSBC( lag = nlags )
	tsbc.fit( Y , X )
	
	lZ = []
	for i in range(nlags + 1):
		tsbc.ref = i
		lZ.append( tsbc.predict(X) )
	
	if plot:
		## Gaussian kernel
		##================
		kde_X = sc.gaussian_kde(X.ravel())
		kde_Y = sc.gaussian_kde(Y.ravel())
		lkde_Z = [ sc.gaussian_kde(Z.ravel()) for Z in lZ ]
		
		## Plot
		##=====
		
		greens = [c for c in plt.cm.Greens( np.linspace(0.3,0.9,nlags+1) ) ]
		
		xmin = min( [T.min() for T in [X,Y] + lZ] )
		xmax = max( [T.max() for T in [X,Y] + lZ] )
		delta = (xmax - xmin) / 10
		xmin -= delta
		xmax += delta
		bins = np.linspace( xmin , xmax , 200 )
		
		nrow,ncol,fs = 1,2,7
		fig = plt.figure( figsize = (fs*ncol,fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.plot( bins , kde_X(bins) , color = "red" , label = r"$X$" )
		for i,kde_Z in enumerate(lkde_Z):
			if i == len(lkde_Z) - 1:
				ax.plot( bins , kde_Z(bins) , color = greens[i] , label = "OTC" )
			else:
				ax.plot( bins , kde_Z(bins) , color = greens[i] )
		ax.plot( bins , kde_Y(bins) , color = "blue" , marker = "" , label = r"$Y$" )
		ax.legend( loc = "upper left" )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		AR.plot_pacf( [X,Y] + lZ , ["red","blue"] + greens , ax = ax , nlags = 7 )
		
		fig.set_tight_layout(True)
		plt.show()
##}}}

def test_schaake_shuffle( plot = True ):##{{{
	X = np.random.normal( size = (10,2) )
	Y = np.random.normal( size = (10,2) )
	
	
	ssr = bct.SchaakeShuffleRef(0,Y)
	Z = ssr.predict(X)
	if plot:
		rX = np.apply_along_axis( sc.rankdata , 0 , X , method = "ordinal" )
		rY = np.apply_along_axis( sc.rankdata , 0 , Y , method = "ordinal" )
		rZ = np.apply_along_axis( sc.rankdata , 0 , Z , method = "ordinal" )
		print(np.hstack( (rX,rY,rZ) ) )
##}}}

def test_metrics():##{{{
	X = np.random.multivariate_normal( mean = np.zeros(2) , cov = np.identity(2) , size = 10000 )
	Y = np.random.multivariate_normal( mean = [5,5]       , cov = np.array( [0.5,-1,-1,2] ).reshape(2,2) , size = 8000 )
	
	bw = bct.bin_width_estimator( [X,Y] )
	muX = bct.SparseHist( X , bw )
	muY = bct.SparseHist( Y , bw )
	
	for metrics in [bcm.chebyshev,bcm.energy,bcm.euclidean,bcm.manhattan,bcm.wasserstein]:
		m = metrics( muX , muY )
		m = metrics( X   , muY )
		m = metrics( muX , Y   )
		m = metrics( X   , Y   )
	
	m = bcm.minkowski( muX , muY , p = 3 )
	m = bcm.minkowski( X   , muY , p = 3 )
	m = bcm.minkowski( muX , Y   , p = 3 )
	m = bcm.minkowski( X   , Y   , p = 3 )
##}}}



def run_all_test( plot = False ):##{{{
	test_qm(          plot = plot )
	test_otc_univ(    plot = plot )
	test_otc_biv(     plot = plot )
	test_ECBC(        plot = plot )
	test_qmrs(        plot = plot )
	test_QDM(         plot = plot )
	test_MBCn(        plot = plot )
	test_MRec(        plot = plot )
#	test_tsbc(        plot = plot )
#	test_metrics()
##}}}


##########
## main ##
##########


if __name__ == "__main__":
	
	print(bc.__version__)
	np.random.seed(42)
	
	
	## Run tests
	##==========
	##{{{
	args_run     = False
	args_verbose = False
	if len(sys.argv) > 1:
		for arg in sys.argv[1:]:
			if arg == "-r" or "arg" == "--run-all-tests":
				args_run = True
			if arg == "-v" or "arg" == "--verbose":
				args_verbose = True
	
	if args_run:
		run_all_test(args_verbose)
	
	##}}}
	
	print( "Done" )

