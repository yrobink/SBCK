
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

#' Manhattan distance
#'
#' Compute Manhattan distance between two dataset or SparseHist X and Y
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
#' d = SBCK::manhattan(  X ,   Y )
#' d = SBCK::manhattan(muX ,   Y )
#' d = SBCK::manhattan(  X , muY )
#' d = SBCK::manhattan(muX , muY )
#'
#' @export
manhattan = function( X , Y )
{
	invisible( SBCK::minkowski(X,Y,1) )
}
