
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

#' Minkowski distance
#'
#' Compute Minkowski distance between two dataset or SparseHist X and Y. If
#' p = 2, it is the Euclidean distance, for p = 1, it is the manhattan distance,
#' if p = Inf, chebyshev distance is called.
#'
#' @param X [matrix or SparseHist] If matrix, dim = ( nrow = n_samples, ncol =
#'        n_features)
#' @param Y [matrix or SparseHist] If matrix, dim = ( nrow = n_samples, ncol =
#'        n_features)
#' @param p [float] power of distance
#'        
#' @return [float] value of distance
#'
#' @examples
#' X = base::cbind( stats::rnorm(2000) , stats::rnorm(2000)  )
#' Y = base::cbind( stats::rnorm(2000,mean=2)  , stats::rnorm(2000) )
#' bw = base::c(0.1,0.1)
#' muX = SBCK::SparseHist( X , bw )
#' muY = SBCK::SparseHist( Y , bw )
#' 
#' ## The four are equals
#' d = SBCK::minkowski(  X ,   Y , p = 3 )
#' d = SBCK::minkowski(muX ,   Y , p = 3 )
#' d = SBCK::minkowski(  X , muY , p = 3 )
#' d = SBCK::minkowski(muX , muY , p = 3 )
#'
#' @export
minkowski = function( X , Y , p = 2 )
{
	if( p == Inf )
		return( SBCK::chebyshev(X,Y) )
	
	mu = SBCK::data_to_hist(X,Y)
	
	muX = mu$muX
	muY = mu$muY
	
	dist = 0
	indx = muY$argwhere( muX$c )
	indy = muX$argwhere( muY$c )
	
	## Common elements of muX in muY
	ii = base::which( indx > 0 )
	if( length(ii) > 0 )
		dist = dist + base::sum( base::abs( muX$p[ii] - muY$p[indx[ii]] )^p )
	
	## Elements of muX not in muY
	ii = base::which(indx == 0)
	if( length(ii) > 0 )
		dist = dist + base::sum( base::abs( muX$p[ ii ] )^p )
	
	## Elements of muY not in muX
	ii = base::which(indy == 0)
	if( length(ii) > 0 )
		dist = dist + base::sum( base::abs( muY$p[ base::which(indy == 0) ] )^p )
	
	dist = dist^(1. / p)
	
	invisible(dist)
}
