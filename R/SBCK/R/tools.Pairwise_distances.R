
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

###############
## Functions ##
###############


#' Pairwise distances
#'
#' Compute the matrix of pairwise distances between a matrix X and a matrix Y
#' 
#' @usage pairwise_distances(X,Y,metric)
#' @param X [matrix] A first matrix (samples in row, features in columns).
#' @param Y [matrix] A second matrix (samples in row, features in columns).
#'        If Y = NULL, then pairwise distances is computed between X and X
#' @param metric [string or callable] The metric used. If metric is a string,
#'        then metric is compiled (so faster). Available string are:
#'        "euclidean", "sqeuclidean" (Square of Euclidean distance),
#'        "logeulidean" (log of the Euclidean distance) and "chebyshev" (max). 
#'        Callable must be a function taking two vectors and returning a double.
#'
#' @return distXY [matrix] Pairwise distances. distXY[i,j] is the distance 
#'         between X[i,] and Y[j,]
#'
#' @examples
#' X = matrix( stats::rnorm(200) , ncol = 100 , nrow = 2 )
#' Y = matrix( stats::rexp(300)  , ncol = 150 , nrow = 2 )
#'
#' distXY = SBCK::pairwise_distances( X , Y ) 
#'
#' @export
pairwise_distances = function( X , Y = NULL , metric = "euclidean" )
{
	if( is.null(Y) )
	{
		if( is.character(metric) )
		{
			return( SBCK::cpp_pairwise_distances_Xstr( X , metric ) )
		}
		if( is.function(metric) )
		{
			return( SBCK::cpp_pairwise_distances_XCall( X , metric ) )
		}
	}
	else
	{
		if( is.character(metric) )
		{
			return( SBCK::cpp_pairwise_distances_XYstr( X , Y , metric ) )
		}
		if( is.function(metric) )
		{
			return( SBCK::cpp_pairwise_distances_XYCall( X , Y , metric ) )
		}
	}
}


