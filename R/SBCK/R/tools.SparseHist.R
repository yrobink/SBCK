
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


#' SparseHist
#'
#' Return the Rcpp Class SparseHistBase initialized
#'
#' @param X  [matrix] Dataset to find the SparseHist
#' @param bin_width   [vector] Width of a bin for each dimension
#' @param bin_origin  [vector] Coordinate of the "0" bin
#'
#' @return [SparseHist] SparseHist class
#'
#' @examples
#' ## Data
#' X = base::matrix( stats::rnorm( n = 10000 ) , nrow = 5000 , ncol = 2 )
#' muX = SparseHist(X)
#'
#' print(muX$p) ## Vector of probabilities
#' print(muX$c) ## Matrix of coordinates of each bins
#' print(muX$argwhere(X)) ## Index of bins of dataset X
#' 
#' @export
SparseHist = function( X , bin_width = NULL , bin_origin = NULL )
{
	## If X is a vector, transform it to matrix
	if( !is.matrix(X) ) X = matrix( X , nrow = length(X) , ncol = 1 )
	
	## bin_width
	if( is.null(bin_width) )
		bin_width = bin_width_estimator(X)
	
	## bin_origin
	if( is.null(bin_origin) )
		bin_origin = base::rep( 0 , length(bin_width) )
	
	muX = SparseHistBase$new( X , bin_width , bin_origin )
	
	invisible( muX )
}
