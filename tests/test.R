
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

base::rm( list = base::ls() )


###############
## Libraries ##
###############

library(R6)
library(devtools)
library(roxygen2)
library(ggplot2)
library(gridExtra)

try(roxygen2::roxygenize("../R/SBCK"))
devtools::load_all("../R/SBCK")
roxygen2::roxygenize("../R/SBCK")
devtools::load_all("../R/SBCK")


###########################
## Useful plot functions ##
###########################

PlotTools = R6::R6Class( "PlotTools" , ##{{{
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	os = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		self$os = self$get_os()
	},
	
	
	#############
	## Methods ##
	#############
	
	get_os = function()
	{
		sysinf = base::Sys.info()
		if( !is.null(sysinf) )
		{
			os = sysinf['sysname']
			if( os == 'Darwin' ) os = "osx"
		}
		else
		{
			## mystery machine
			os = .Platform$OS.type
			if( base::grepl( "^darwin"   , R.version$os ) ) os = "osx"
			if( base::grepl( "linux-gnu" , R.version$os ) ) os = "linux"
		}
		invisible(tolower(os))
	},
	
	new_screen = function()
	{
		if( self$os == "osx" )
		{
			grDevices::quartz()
		}
		if( self$os == "linux" )
		{
			grDevices::X11()
		}
	},
	
	wait = function()
	{
		while( base::names(grDevices::dev.cur()) !='null device' ) base::Sys.sleep(1)
	}
	
	)
)
##}}}

plt = PlotTools$new()


###############
## Functions ##
###############

## Some functions
##===============

## Test tools
##===========

test_pairwise_distances = function( show = FALSE )##{{{
{
	X = matrix( stats::rnorm(200) , nrow = 100 , ncol = 2 )
	Y = matrix( stats::rexp(200)  , nrow = 100 , ncol = 2 )
	
	distXY  = SBCK::pairwise_distances( X , Y ) 
	distX   = SBCK::pairwise_distances(X)
	ldistXY = SBCK::pairwise_distances( X , Y , "sqeuclidean" )
	ndistXY = SBCK::pairwise_distances( X , Y , function( x , y ) { return(-1) } ) ## Stupid distance, just to check if == -1 for all values
}
##}}}

test_rv_histogram = function( show = FALSE )##{{{
{
	X = stats::rnorm( n = 10000 )
	
	rvX = SBCK::rv_histogram$new(X)
	
	## RVS
	X2 = rvX$rvs(10000)
	
	## CDF and survival function (1-CDF)
	x   = base::seq( rvX$min , rvX$max , length = 1000 )
	cdf = rvX$cdf(x)
	sf  = rvX$sf(x)
	
	## cdf inverse and sf inverse
	p = base::seq( 0.05 , 0.95 , length = 1000 )
	icdf = rvX$icdf(p)
	isf  = rvX$isf(p)
	
	if( show )
	{
		plt$new_screen()
		graphics::par( mfrow = base::c( 2 , 3 ) )
		
		graphics::plot( x , cdf  , col = "red" , main = "CDF" )
		graphics::plot( x , sf   , col = "red" , main = "sf" )
		graphics::hist( X , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , freq = FALSE )
		graphics::lines( density(X)  , col = "red" )
		graphics::plot( p , icdf , col = "red" , main = "iCDF" )
		graphics::plot( p , isf  , col = "red" , main = "isf" )
		graphics::hist( X2 , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , freq = FALSE )
		graphics::lines( density(X2)  , col = "red" )
		
	}
}
##}}}

test_SparseHist = function( show = FALSE )##{{{
{
	X = base::cbind( stats::rnorm(100000) , stats::rexp(100000) )
	muX = SBCK::SparseHist(X)
	
	if(show)
	{
		plt$new_screen()
		g = ggplot2::ggplot( data.frame( x = muX$c[,1] , y = muX$c[,2] ) , ggplot2::aes( x = x , y = y ) )
		g = g + ggplot2::geom_point( ggplot2::aes( color = muX$p ) )
		g = g + ggplot2::scale_color_viridis_c()
		g = g + ggplot2::ggtitle( base::paste( "SparseHist of normal x exp," , muX$n_samples , "bins found," , muX$n_features , "features." ) )
		graphics::plot(g)
	}
}
##}}}

test_OT = function( show = FALSE )##{{{
{
	X = stats::rnorm(2000)
	Y = stats::rnorm(2000 , mean = 5 )
	bw = base::c(0.1)
	muX = SBCK::SparseHist( X , bw )
	muY = SBCK::SparseHist( Y , bw )
	distXY = SBCK::pairwise_distances( muX$c , muY$c )
	
	ot = OTNetworkSimplex$new()
	ot$fit( muX , muY )
	
	if(show)
	{
		cat( "Test OT\n" )
		print( base::paste( "sum of plan =" , sum(ot$plan) ) )
		print( base::paste( "success =" , ot$success ) )
		print( base::paste( "cost =" , sqrt(sum(ot$plan * ot$C)) ) )
	}
}
##}}}

## Test metrics
##=============

test_metrics = function( show = FALSE )##{{{
{
	## Data
	X = base::cbind( stats::rnorm(2000) , stats::rnorm(2000)  )
	Y = base::cbind( stats::rnorm(2000,mean=10)  , stats::rnorm(2000) )
	
	bw = base::c(0.1,0.1)
	muX = SBCK::SparseHist( X , bw )
	muY = SBCK::SparseHist( Y , bw )
	
	lmetric = base::c( SBCK::wasserstein , SBCK::energy , SBCK::chebyshev , SBCK::euclidean , SBCK::manhattan )
	lname   = base::c( "wasserstein" , "energy" , "chebyshev" , "euclidean" , "manhattan" )
	
	i = 1
	for( metric in lmetric )
	{
		d0 = metric( X   ,   Y )
		d1 = metric( X   , muY )
		d2 = metric( muX ,   Y )
		d3 = metric( muX , muY )
		
		if(show)
		{
			print(lname[i])
			print(d0)
			print(d1)
			print(d2)
			print(d3)
		}
		i = i + 1
	}
}
##}}}


## Test BC method
##===============

test_qm = function( show = FALSE )##{{{
{
	X0 = stats::rnorm( n = 2000 , mean = 0 , sd = 1   )
	X1 = stats::rnorm( n = 2001 , mean = 5 , sd = 2   )
	Y0 = stats::rnorm( n = 2002 , mean = 2 , sd = 0.5 )
	
	qm = SBCK::QM$new()
	qm$fit( Y0 , X0 )
	Z0 = qm$predict(X0)
	
	cdft = SBCK::CDFt$new()
	cdft$fit( Y0 , X0 , X1 )
	Z1 = cdft$predict( X1 )
	
	if(show)
	{
		plt$new_screen()
		
		
		g0 = ggplot2::ggplot() + ggplot2::aes( x = x )
		g0 = g0 + ggplot2::geom_histogram( binwidth = 0.1 , aes( y = ..density.. ) , data = data.frame(x = X0) , color = "red"  , alpha = 0.5 )
		g0 = g0 + ggplot2::geom_histogram( binwidth = 0.1 , aes( y = ..density.. ) , data = data.frame(x = Y0) , color = "blue" , alpha = 0.5 )
		g0 = g0 + ggplot2::ggtitle( "Calibration" )
		
		g1 = ggplot2::ggplot() + ggplot2::aes( x = x )
		g1 = g1 + ggplot2::geom_histogram( binwidth = 0.1 , aes( y = ..density.. ) , data = data.frame(x = X0) , color = "red"   , alpha = 0.5 )
		g1 = g1 + ggplot2::geom_histogram( binwidth = 0.1 , aes( y = ..density.. ) , data = data.frame(x = Z0) , color = "green" , alpha = 0.5 )
		g1 = g1 + ggplot2::ggtitle( "QM" )
		
		g2 = ggplot2::ggplot() + ggplot2::aes( x = x )
		g2 = g2 + ggplot2::geom_histogram( binwidth = 0.1 , aes( y = ..density.. ) , data = data.frame(x = X0) , color = "red"  , alpha = 0.5 )
		g2 = g2 + ggplot2::geom_histogram( binwidth = 0.1 , aes( y = ..density.. ) , data = data.frame(x = Y0) , color = "blue" , alpha = 0.5 )
		g2 = g2 + ggplot2::ggtitle( "Calibration" )
		
		g3 = ggplot2::ggplot() + ggplot2::aes( x = x )
		g3 = g3 + ggplot2::geom_histogram( binwidth = 0.1 , aes( y = ..density.. ) , data = data.frame(x = X1) , color = "red"   , alpha = 0.5 )
		g3 = g3 + ggplot2::geom_histogram( binwidth = 0.1 , aes( y = ..density.. ) , data = data.frame(x = Z1) , color = "green" , alpha = 0.5 )
		g3 = g3 + ggplot2::ggtitle( "CDFt" )
		
		g = gridExtra::grid.arrange( g0 , g1 , g2 , g3 , ncol = 2 , nrow = 2 )
		
		graphics::plot(g)
	}
}
##}}}

test_qmrs = function( show = FALSE )##{{{
{
	X0 = base::cbind( stats::rnorm(2000) , stats::rexp(2000)  )
	Y0 = base::cbind( stats::rexp(2000)  , stats::rnorm(2000) )
	X1 = base::cbind( stats::rnorm(2000 , mean = 5 ) , stats::rexp(2000)  )
	
	qmrs = SBCK::QMrs$new()
	qmrs$fit( Y0 , X0 )
	Z0 = qmrs$predict(X0)
	
	r2d2 = SBCK::R2D2$new()
	r2d2$fit( Y0 , X0 , X1 )
	Z1 = r2d2$predict(X1)
	
	if(show)
	{
		plt$new_screen()
		
		g0 = ggplot2::ggplot( data.frame( x = X0[,1] , y = X0[,2] ) , ggplot2::aes( x = x , y = y ) )
		g0 = g0 + ggplot2::geom_point( color = "red" )
		g0 = g0 + ggplot2::ggtitle( "X0" )
		
		g1 = ggplot2::ggplot( data.frame( x = Y0[,1] , y = Y0[,2] ) , ggplot2::aes( x = x , y = y ) )
		g1 = g1 + ggplot2::geom_point( color = "blue" )
		g1 = g1 + ggplot2::ggtitle( "Y0" )
		
		g2 = ggplot2::ggplot( data.frame( x = Z0[,1,1] , y = Z0[,2,1] ) , ggplot2::aes( x = x , y = y ) )
		g2 = g2 + ggplot2::geom_point( color = "green" )
		g2 = g2 + ggplot2::ggtitle( "Z0" )
		
		g3 = ggplot2::ggplot( data.frame( x = X1[,1] , y = X1[,2] ) , ggplot2::aes( x = x , y = y ) )
		g3 = g3 + ggplot2::geom_point( color = "red" )
		g3 = g3 + ggplot2::ggtitle( "X1" )
		
		g5 = ggplot2::ggplot( data.frame( x = Z1[,1,1] , y = Z1[,2,1] ) , ggplot2::aes( x = x , y = y ) )
		g5 = g5 + ggplot2::geom_point( color = "green" )
		g5 = g5 + ggplot2::ggtitle( "Z1" )
		
		
		g = gridExtra::grid.arrange( g0 , g1 , g2 , g3 , ggplot2::ggplot() , g5 , ncol = 3 , nrow = 2 )
		graphics::plot(g)
	}
}
##}}}

test_otc = function( show = FALSE )##{{{
{
	X0 = base::cbind( stats::rnorm(2000) , stats::rexp(2000)  )
	Y0 = base::cbind( stats::rexp(2000)  , stats::rnorm(2000) )
	X1 = base::cbind( stats::rnorm(2000 , mean = 5 ) , stats::rexp(2000)  )
	
	otc = SBCK::OTC$new()
	otc$fit( Y0 , X0 )
	Z0 = otc$predict(X0)
	
	dotc = SBCK::dOTC$new( cov_factor = "cholesky" )
	dotc$fit( Y0 , X0 , X1 )
	Z1 = dotc$predict(X1)
	
	if( show )
	{
		plt$new_screen()
		
		g0 = ggplot2::ggplot( data.frame( x = X0[,1] , y = X0[,2] ) , ggplot2::aes( x = x , y = y ) )
		g0 = g0 + ggplot2::geom_point( color = "red" )
		g0 = g0 + ggplot2::ggtitle( "X0" )
		
		g1 = ggplot2::ggplot( data.frame( x = Y0[,1] , y = Y0[,2] ) , ggplot2::aes( x = x , y = y ) )
		g1 = g1 + ggplot2::geom_point( color = "blue" )
		g1 = g1 + ggplot2::ggtitle( "Y0" )
		
		g2 = ggplot2::ggplot( data.frame( x = Z0[,1] , y = Z0[,2] ) , ggplot2::aes( x = x , y = y ) )
		g2 = g2 + ggplot2::geom_point( color = "green" )
		g2 = g2 + ggplot2::ggtitle( "Z0" )
		
		g3 = ggplot2::ggplot( data.frame( x = X1[,1] , y = X1[,2] ) , ggplot2::aes( x = x , y = y ) )
		g3 = g3 + ggplot2::geom_point( color = "red" )
		g3 = g3 + ggplot2::ggtitle( "X1" )
		
		g5 = ggplot2::ggplot( data.frame( x = Z1[,1] , y = Z1[,2] ) , ggplot2::aes( x = x , y = y ) )
		g5 = g5 + ggplot2::geom_point( color = "green" )
		g5 = g5 + ggplot2::ggtitle( "Z1" )
		
		g = gridExtra::grid.arrange( g0 , g1 , g2 , g3 , ggplot2::ggplot() , g5 , ncol = 3 , nrow = 2 )
		graphics::plot(g)
	}
}
##}}}


## All in one
##===========

run_all_tests = function( show = FALSE )##{{{
{
	## Tools tests
	test_pairwise_distances(show)
	test_rv_histogram(show)
	test_SparseHist(show)
	test_OT(show)
	
	## Metrics tests
	test_metrics(show)
	
	## BC tests
	test_qm(show)
	test_qmrs(show)
	test_otc(show)
}
##}}}



##########
## main ##
##########

## Read command line arguments and run (or not) tests
##================================================{{{

args = commandArgs( trailingOnly = TRUE )
args_verbose = FALSE
args_run     = FALSE
if( length(args) > 0 )
{
	for( a in args )
	{
		if( a == "-r" || a == "--run-all-tests" )
			args_run = TRUE
		if( a == "-v" || a == "--verbose" )
			args_verbose = TRUE
	}
}

if( args_run )
	run_all_tests(args_verbose)

##}}}

plt$wait()
print("Done")

