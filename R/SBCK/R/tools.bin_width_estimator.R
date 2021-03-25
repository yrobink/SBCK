
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

#' bin_width_estimator method
#'
#' Lenght of cell to compute an histogram
#'
#' @param X [matrix] A matrix containing data, nrow = n_samples,
#'        ncol = n_features
#' @param method [string] Method to estimate bin_width, values are "auto", "FD"
#'        (Friedman Draconis, robust over outliners) or "Sturges". If "auto" is
#'        used and if nrow(X) < 1000, "Sturges" is used, else "FD" is used.
#'        
#' @return [vector] Lenght of bins
#'
#' @importFrom stats quantile
#' @examples
#' X = base::cbind( stats::rnorm( n = 2000 ) , stats::rexp(2000) )
#' ## Friedman Draconis is used
#' binw_width = SBCK::bin_width_estimator( X , method = "auto" ) 
#' X = stats::rnorm( n = 500 )
#' ## Sturges is used
#' binw_width = SBCK::bin_width_estimator( X , method = "auto" ) 
#' 
#' @export
bin_width_estimator = function( X , method = "auto" ) 
{
	if( "list" %in% class(X) )
	{
		bw = matrix( NA , nrow = length(X) , ncol = base::ncol(X[[1]]) )
		for( i in 1:length(X) )
			bw[i,] = SBCK::bin_width_estimator( X[[i]] , method )
		return( base::apply( bw , 2 , base::min ) )
	}
	
	if( !is.matrix(X) )
		X = matrix( X , nrow = length(X) , ncol = 1 )
	
	n_samples  = dim(X)[1]
	n_features = dim(X)[2]
	
	## Find method to use
	if( method == "auto" && n_samples < 1000 )
	{
		method = "Sturges"
	}
	
	## Find bin_width
	bin_width = rep( 0 , n_features )
	if( method == "Sturges" )
	{
		nh = log2(n_samples) + 1
		bin_width = rep( 1. / nh , n_features )
	}
	else ## FD (Freedman Diaconis) method, robust over outliers
	{
		pow = n_samples^(1./3.)
		q    = base::apply( X , 2 , quantile , probs = base::c( 0.25 , 0.75 ) )
		bin_width = 2 * ( q[2,] - q[1,] ) / pow
	}
	
	return(as.vector(bin_width))
}

