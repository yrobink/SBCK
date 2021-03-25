
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

base::rm( list = base::ls() )


###############
## Libraries ##
###############

library(R6)
library(devtools)
library(roxygen2)
library(ggplot2)
library(gridExtra)
library(ROOPSD)
library(pmetric)

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
	distXY = pmetric::pairwise_distances( muX$c , muY$c )
	
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

test_shuffle = function( show = FALSE )##{{{
{
	X0 = matrix( stats::runif(30) , nrow = 10 , ncol = 3 )
	Y0 = matrix( stats::runif(30) , nrow = 10 , ncol = 3 )
	
	ssr = SBCK::SchaakeShuffleRef$new( ref = 1 )
	ssr$fit(Y0)
	Z0 = ssr$predict(X0)
	
	rank_X0 = base::apply( X0 , 2 , base::rank )
	rank_Y0 = base::apply( Y0 , 2 , base::rank )
	rank_Z0 = base::apply( Z0 , 2 , base::rank )
	
	if(show)
	{
		cat( "Test shuffle\n" )
		print( base::cbind( rank_X0 , rank_Y0 , rank_Z0 ) )
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
	XY = SBCK::dataset_gaussian_exp_mixture_1d(2000)
	X0 = XY$X0
	X1 = XY$X1
	Y0 = XY$Y0
	
	qm = SBCK::QM$new()
	qm$fit( Y0 , X0 )
	Z0 = qm$predict(X0)
	
	cdft = SBCK::CDFt$new()
	cdft$fit( Y0 , X0 , X1 )
	Z1 = cdft$predict(X1)
	
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

test_ECBC = function( show = FALSE )##{{{
{
	XY = SBCK::dataset_bimodal_reverse_2d(2000)
	X0 = XY$X0
	X1 = XY$X1
	Y0 = XY$Y0
	
	ecbc = SBCK::ECBC$new()
	ecbc$fit( Y0 , X0 , X1 )
	Z1 = ecbc$predict(X1)
	Z0 = ecbc$predict(X1,X0)$Z0
	
	if(show)
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

test_qmrs = function( show = FALSE )##{{{
{
	XY = SBCK::dataset_bimodal_reverse_2d(2000)
	X0 = XY$X0
	X1 = XY$X1
	Y0 = XY$Y0
	
	r2d2 = SBCK::R2D2$new()
	r2d2$fit( Y0 , X0 , X1 )
	Z1 = r2d2$predict(X1)
	Z0 = r2d2$predict(X1,X0)$Z0
	
	if(show)
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

test_otc_univ = function( show = FALSE )##{{{
{
	XY = SBCK::dataset_gaussian_exp_mixture_1d(2000)
	X0 = XY$X0
	X1 = XY$X1
	Y0 = XY$Y0
	
	dotc = SBCK::dOTC$new()
	dotc$fit( Y0 , X0 , X1 )
	Z1 = dotc$predict(X1)
	Z0 = dotc$predict( X1 , X0 )$Z0
	
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

test_otc = function( show = FALSE )##{{{
{
	XY = SBCK::dataset_bimodal_reverse_2d(2000)
	X0 = XY$X0
	X1 = XY$X1
	Y0 = XY$Y0
	
	dotc = SBCK::dOTC$new()
	dotc$fit( Y0 , X0 , X1 )
	Z1 = dotc$predict(X1)
	Z0 = dotc$predict(X1,X0)$Z0
	
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

test_MRec = function( show = FALSE )##{{{
{
	XY = SBCK::dataset_bimodal_reverse_2d(2000)
	X0 = XY$X0
	X1 = XY$X1
	Y0 = XY$Y0
	
	mrec = SBCK::MRec$new()
	mrec$fit( Y0 , X0 , X1 )
	Z1 = mrec$predict(X1)
	Z0 = mrec$predict(X1,X0)$Z0
	
	
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

test_QDM = function( show = FALSE )##{{{
{
	XY = SBCK::dataset_bimodal_reverse_2d(2000)
	X0 = XY$X0
	X1 = XY$X1
	Y0 = XY$Y0
	
	qdm = SBCK::QDM$new()
	qdm$fit( Y0 , X0 , X1 )
	Z1 = qdm$predict(X1)
	Z0 = qdm$predict(X1,X0)$Z0
	
	
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

test_MBCn = function( show = FALSE )##{{{
{
	XY = SBCK::dataset_bimodal_reverse_2d(2000)
	X0 = XY$X0
	X1 = XY$X1
	Y0 = XY$Y0
	
	mbcn = SBCK::MBCn$new()
	mbcn$fit( Y0 , X0 , X1 )
	Z1 = mbcn$predict(X1)
	Z0 = mbcn$predict(X1,X0)$Z0
	
	
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
		
		
		nit = mbcn$iter_slope$nit - 1
		g4 = ggplot2::ggplot() + ggplot2::geom_line( ggplot2::aes( x = 1:nit , y = mbcn$iter_slope$criteria[1:nit] ) )
		g4 = g4 + ggplot2::ggtitle( "MBCn convergence" ) + ggplot2::labs( x = "Iterations" , y = "Wasserstein" )
		
		g5 = ggplot2::ggplot( data.frame( x = Z1[,1] , y = Z1[,2] ) , ggplot2::aes( x = x , y = y ) )
		g5 = g5 + ggplot2::geom_point( color = "green" )
		g5 = g5 + ggplot2::ggtitle( "Z1" )
		
		g = gridExtra::grid.arrange( g0 , g1 , g2 , g3 , g4 , g5 , ncol = 3 , nrow = 2 )
		graphics::plot(g)
	}
}
##}}}


## All in one
##===========

run_all_tests = function( show = FALSE )##{{{
{
	## Tools tests
	test_SparseHist(show)
	test_OT(show)
	test_shuffle(show)
	
	## Metrics tests
	test_metrics(show)
	
	## BC tests
	test_qm(show)
	test_ECBC(show)
	test_qmrs(show)
	test_otc_univ(show)
	test_otc(show)
	test_MRec(show)
	test_QDM(show)
	test_MBCn(show)
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

