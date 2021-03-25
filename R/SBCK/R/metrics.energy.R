
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

#' Energy distance
#'
#' Compute Energy distance between two dataset or SparseHist X and Y
#'
#' @param X [matrix or SparseHist] If matrix, dim = ( nrow = n_samples, ncol =
#'        n_features)
#' @param Y [matrix or SparseHist] If matrix, dim = ( nrow = n_samples, ncol =
#'        n_features)
#' @param p [float] power of energy distance, default is 2.
#' @param metric [str or function] metric for pairwise distance, default is
#'        "euclidean", see SBCK::pairwise_distances
#'        
#' @return [float] value of distance
#'
#' @examples
#' X = base::cbind( stats::rnorm(2000) , stats::rnorm(2000)  )
#' Y = base::cbind( stats::rnorm(2000,mean=10)  , stats::rnorm(2000) )
#' bw = base::c(0.1,0.1)
#' muX = SBCK::SparseHist( X , bw )
#' muY = SBCK::SparseHist( Y , bw )
#' 
#' ## The four are equals
#' w2 = SBCK::energy(X,Y)
#' w2 = SBCK::energy(muX,Y)
#' w2 = SBCK::energy(X,muY)
#' w2 = SBCK::energy(muX,muY)
#'
#' @export
energy = function( X , Y , p = 2 , metric = "euclidean" )
{
	mu = SBCK::data_to_hist(X,Y)
	muX = mu$muX
	muY = mu$muY
	
	pX = matrix( muX$p , nrow = muX$n_samples , ncol = 1 )
	pY = matrix( muY$p , nrow = muY$n_samples , ncol = 1 )
	
	XY = base::sum( SBCK::pairwise_distances( muX$c , muY$c , metric = metric )^p * ( pX %*% base::t(pY) ) )
	XX = base::sum( SBCK::pairwise_distances( muX$c         , metric = metric )^p * ( pX %*% base::t(pX) ) )
	YY = base::sum( SBCK::pairwise_distances( muY$c         , metric = metric )^p * ( pY %*% base::t(pY) ) )
	
	cost = (2 * XY - XX - YY)^( 1. / p )
	invisible(cost)
}
