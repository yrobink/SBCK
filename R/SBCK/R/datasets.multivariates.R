
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

# dataset_gaussian_exp_2d ##{{{

#' dataset_gaussian_exp_2d
#'
#' Generate a testing dataset such that the biased dataset is a distribution
#' of the the form Normal x Exp and the reference of the the form Exp x Normal.
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection) 
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_gaussian_exp_2d = function(n_samples)
{
	X0 = base::cbind( stats::rnorm(n_samples)             , stats::rexp( n_samples)  )
	Y0 = base::cbind( stats::rexp( n_samples)             , stats::rnorm(n_samples) )
	X1 = base::cbind( stats::rnorm(n_samples , mean = 5 ) , stats::rexp( n_samples)  )
	
	return( list( X0 = X0 , X1 = X1 , Y0 = Y0 ) )
}
##}}}

# dataset_gaussian_L_2d ##{{{

#' dataset_gaussian_L_2d
#'
#' Generate a testing dataset such that the biased dataset is a normal
#' distribution and reference a mixture a normal with a form in "L"
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection)
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_gaussian_L_2d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_gaussian_L_2d = function( n_samples )
{
	## Construction of X0 (biased period 0), X1 (biased period 1) and Y0 (reference period 0)
	size0  = as.integer(n_samples/2)
	size1  = n_samples - as.integer(n_samples/4)
	
	## Just a gaussian for X0
	X0 = ROOPSD::rmultivariate_normal( n = n_samples , mean = base::c(0.,0.) , cov = base::diag(2) )
	
	## A lightly complex gaussian for X1
	X1 = ROOPSD::rmultivariate_normal( n = n_samples , mean = base::c(1.,2.) , cov = matrix( base::c(2.,0,0,0.5) , nrow = 2 , ncol = 2 ) )
	
	## A very complex law for Y0
	Y0 = matrix( NA , nrow = n_samples , ncol = 2 )
	Y0[1:size0,]         = ROOPSD::rmultivariate_normal( n = size0             , mean = base::c(7.,7.)   , cov = matrix( base::c(2,0,0,0.5)   , nrow = 2 , ncol = 2 ) )
	Y0[size0:size1,]     = ROOPSD::rmultivariate_normal( n = size1 - size0 + 1 , mean = base::c(5.,9.)   , cov = matrix( base::c(0.5,0,0,2)   , nrow = 2 , ncol = 2 ) )
	Y0[size1:n_samples,] = ROOPSD::rmultivariate_normal( n = n_samples - size1 + 1 , mean = base::c(5.,12.5) , cov = matrix( base::c(0.2,0,0,0.2) , nrow = 2 , ncol = 2 ) )
	meanY0 = base::apply( Y0 , 2 , base::mean )
	meanX0 = base::apply( X0 , 2 , base::mean )
	diff = meanY0 - meanX0
	Y0 = base::t(base::apply( Y0 , 1 , function(x) { return(x - diff) } ))
	
	return( list( X0 = X0 , X1 = X1 , Y0 = Y0 ) )
}
##}}}

# dataset_gaussian_2d ##{{{

#' dataset_gaussian_2d
#'
#' Generate a testing dataset from random bivariate Gaussian distribution
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection)
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_gaussian_2d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_gaussian_2d = function(n_samples)
{
	CX0 = base::matrix( base::c(2.33,-0.38,-0.38,0.48) , nrow = 2 , ncol = 2 )
	X0  = ROOPSD::rmultivariate_normal( n_samples , mean = base::rep(0,2) , cov = CX0 )
	CX1 = base::matrix( base::c(0.1,-0.8,-0.8,3.37) , nrow = 2 , ncol = 2 )
	X1  = ROOPSD::rmultivariate_normal( n_samples , mean = base::rep(2,2) , cov = CX1 )
	CY0 = base::matrix( base::c(0.32,-0.18,-0.18,2.40) , nrow = 2 , ncol = 2 )
	Y0  = ROOPSD::rmultivariate_normal( n_samples , mean = base::rep(3,2) , cov = CY0 )
	
	return( list( X0 = X0 , X1 = X1 , Y0 = Y0 ) )
}
##}}}

# dataset_bimodal_reverse_2d ##{{{

#' dataset_bimodal_reverse_2d
#'
#' Generate a testing dataset from bimodale random bivariate Gaussian distribution
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection)
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_bimodal_reverse_2d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_bimodal_reverse_2d = function(n_samples)
{
	draw   = list( u = as.integer(n_samples/2) )
	draw[["l"]] = n_samples - draw$u
	lmY0   = list( u = base::c(5,-3)  , l = base::c(-3,3) )
	lcovY0 = list( u = 0.9 * diag(2)  , l = diag(2)       )
	lmX0   = list( u = base::c(0,0)   , l = base::c(2,2)  )
	lcovX0 = list( u = diag(2)        , l = 0.5 * diag(2) )
	lmX1   = list( u = base::c(-1,-1) , l = base::c(5,5)  )
	lcovX1 = list( u = 2 * diag(2)    , l = 0.1 * diag(2) )
	
	Y0 = NULL
	X0 = NULL
	X1 = NULL
	for( idx in base::c("u","l") )
	{
		Y0 = base::rbind( Y0 , ROOPSD::rmultivariate_normal( draw[[idx]] , lmY0[[idx]] , lcovY0[[idx]] ) )
		X0 = base::rbind( X0 , ROOPSD::rmultivariate_normal( draw[[idx]] , lmX0[[idx]] , lcovX0[[idx]] ) )
		X1 = base::rbind( X1 , ROOPSD::rmultivariate_normal( draw[[idx]] , lmX1[[idx]] , lcovX1[[idx]] ) )
	}
	
	return( list( X0 = X0 , X1 = X1 , Y0 = Y0 ) )
}

##}}}

