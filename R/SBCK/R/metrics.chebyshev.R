
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

#' Chebyshev distance
#'
#' Compute Chebyshev distance between two dataset or SparseHist X and Y
#'
#' @param X [matrix or SparseHist] If matrix, dim = ( nrow = n_samples, ncol =
#'        n_features)
#' @param Y [matrix or SparseHist] If matrix, dim = ( nrow = n_samples, ncol =
#'        n_features)
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
#' d = SBCK::chebyshev(  X ,   Y )
#' d = SBCK::chebyshev(muX ,   Y )
#' d = SBCK::chebyshev(  X , muY )
#' d = SBCK::chebyshev(muX , muY )
#'
#' @export
chebyshev = function( X , Y )
{
	mu = SBCK::data_to_hist(X,Y)
	
	muX = mu$muX
	muY = mu$muY
	
	dist = 0
	indx = muY$argwhere( muX$c )
	indy = muX$argwhere( muY$c )
	
	## Common elements of muX in muY
	ii = base::which( indx > 0 )
	if( length(ii) > 0 )
		dist = base::max( dist , base::max( base::abs( muX$p[ii] - muY$p[indx[ii]] ) ) )
	
	## Elements of muX not in muY
	ii = base::which(indx == 0)
	if( length(ii) > 0 )
		dist = base::max( dist , base::max( base::abs( muX$p[ ii ] ) ) )
	
	## Elements of muY not in muX
	ii = base::which(indy == 0)
	if( length(ii) > 0 )
		dist = base::max( dist , base::max( base::abs( muY$p[ ii ] ) ) )
	
	invisible(dist)
}
