
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

## dataset_gaussian_exp_mixture_1d ##{{{

#' dataset_gaussian_exp_mixture_1d
#'
#' Generate a univariate testing dataset from a mixture of gaussian and
#' exponential distribution
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection)
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_gaussian_exp_mixture_1d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_gaussian_exp_mixture_1d = function( n_samples )
{
	dsize = as.integer(0.75 * n_samples)
	
	X0 = matrix( NA , nrow = n_samples , ncol = 1 )
	X0[,1] = stats::rnorm( n = n_samples , mean = 7 , sd = 1 )
	
	X1 = matrix( NA , nrow = n_samples , ncol = 1 )
	X1[1:(n_samples-dsize),1]         = stats::rnorm( n = n_samples - dsize , mean = 5 , sd = 1 )
	X1[(n_samples-dsize):n_samples,1] = stats::rnorm( n = dsize + 1 , mean = 9 , sd = 1 )
	
	Y0 = matrix( NA , nrow = n_samples , ncol = 1 )
	Y0[1:dsize,1]         = stats::rexp( n = dsize , rate = 1 )
	Y0[dsize:n_samples,1] = stats::rnorm( n = n_samples - dsize + 1 , mean = 10 , sd = 1 )
	
	return( list( Y0 = Y0 , X0 = X0 , X1 = X1 ) )
}
##}}}

## dataset_gaussian_VS_exp_1d ##{{{

#' dataset_gaussian_VS_exp_1d
#'
#' Generate a univariate testing dataset such that biased data follow an
#' exponential law whereas reference follow a normal distribution
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection)
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_gaussian_VS_exp_1d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_gaussian_VS_exp_1d = function( n_samples )
{
	Y0 = matrix( stats::rnorm( n = n_samples , mean = 5 , sd = 1 ) , nrow = n_samples , ncol = 1 )
	X0 = matrix( stats::rexp(  n = n_samples , rate = 1 ) , nrow = n_samples , ncol = 1 )
	X1 = matrix( stats::rexp(  n = n_samples , rate = 2 ) , nrow = n_samples , ncol = 1 )
	
	return( list( Y0 = Y0 , X0 = X0 , X1 = X1 ) )
}
##}}}


