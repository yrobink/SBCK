
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

#' wasserstein distance
#'
#' Compute wasserstein distance between two dataset or SparseHist X and Y
#'
#' @param X [matrix or SparseHist] If matrix, dim = ( nrow = n_samples, ncol =
#'        n_features)
#' @param Y [matrix or SparseHist] If matrix, dim = ( nrow = n_samples, ncol =
#'        n_features)
#' @param p [float]
#'        Power of the metric (default = 2)
#' @param ot [Optimal transport solver]
#'        
#' @return [float] value of distance
#'
#' @references Wasserstein, L. N. (1969). Markov processes over denumerable
#'             products of spaces describing large systems of automata. Problems
#'             of Information Transmission, 5(3), 47-52.
#'
#' @examples
#' X = base::cbind( stats::rnorm(2000) , stats::rnorm(2000)  )
#' Y = base::cbind( stats::rnorm(2000,mean=10)  , stats::rnorm(2000) )
#' bw = base::c(0.1,0.1)
#' muX = SBCK::SparseHist( X , bw )
#' muY = SBCK::SparseHist( Y , bw )
#' 
#' ## The four are equals
#' w2 = SBCK::wasserstein(X,Y)
#' w2 = SBCK::wasserstein(muX,Y)
#' w2 = SBCK::wasserstein(X,muY)
#' w2 = SBCK::wasserstein(muX,muY)
#'
#' @export
wasserstein = function( X , Y , p = 2 , ot = SBCK::OTNetworkSimplex$new() )
{
	mu = SBCK::data_to_hist(X,Y)
	ot$p = p
	ot$fit( mu$muX , mu$muY )
	
	cost = base::sum(ot$plan * ot$C)^(1. / ot$p)
	if( !ot$success )
	{
		cost = NaN
	}
	if( !( abs( base::sum(ot$plan) - 1 ) < 1e-6 ) )
	{
		cost = NaN
	}
	
	invisible(cost)
}
