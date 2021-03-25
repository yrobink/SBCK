
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

#' data_to_hist
#'
#' Just a function to transform two datasets into SparseHist, if X or Y (or the
#' both) are already a SparseHist, update just the second
#'
#' @param X [matrix or SparseHist]
#' @param Y [matrix or SparseHist]
#'        
#' @return [list(muX,muY)] a list with the two SparseHist
#'
#' @examples
#' X = base::cbind( stats::rnorm(2000) , stats::rexp(2000)  )
#' Y = base::cbind( stats::rexp(2000)  , stats::rnorm(2000) )
#' 
#' bw = base::c(0.1,0.1)
#' muX = SBCK::SparseHist( X , bw )
#' muY = SBCK::SparseHist( Y , bw )
#' 
#' ## The four give the same result
#' SBCK::data_to_hist( X   , Y )
#' SBCK::data_to_hist( muX , Y )
#' SBCK::data_to_hist( X   , muY )
#' SBCK::data_to_hist( muX , muY )
#'
#' @export
data_to_hist = function( X , Y )
{
	is_hist = function(Z) { return( ( "Rcpp_SparseHistBase" %in% class(Z) ) || ("OTHist" %in% class(Z))  ) }
	X_is_hist = is_hist(X)
	Y_is_hist = is_hist(Y)
	
	if( X_is_hist && Y_is_hist )
	{
		return( list( muX = X , muY = Y ) )
	}
	if( X_is_hist && !Y_is_hist )
	{
		muY = SBCK::SparseHist( Y , X$bin_width , X$bin_width )
		return( list( muX = X , muY = muY ) )
	}
	if( !X_is_hist && Y_is_hist )
	{
		muX = SBCK::SparseHist( X , Y$bin_width , Y$bin_width )
		return( list( muX = muX , muY = Y ) )
	}
	
	bw = SBCK::bin_width_estimator( list(X,Y) )
	muX = SBCK::SparseHist( X , bw )
	muY = SBCK::SparseHist( Y , bw )
	
	return( list( muX = muX , muY = muY ) )
}
