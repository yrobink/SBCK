
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


#' Optimal Transport Histogram
#'
#' @description
#' Histogram
#'
#' @details
#' Just a generic class which contains two arguments, p (probability) and c
#' (center of bins)
#'
#' @examples
#' ## Build a random discrete probability distribution
#' p = stats::rnorm(100)
#' p = p / base::sum(p)
#' c = base::seq( -1 , 1 , length = 100 )
#' mu = OTHist$new( p , c )
#'
#' @export
OTHist = R6::R6Class( "OTHist" ,
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field p [vector] Vector of probability
	p          = NULL,
	#' @field c [matrix] Vector of center of bins, with nrow = n_samples and ncol = n_features
	c          = NULL,
	#' @field bin_width [vector or NULL] A vector of lengths of the cells
	#'        discretizing R^{numbers of variables}. If NULL, it is estimating
	#'        during the fit
	bin_width  = NULL,
	#' @field bin_origin [vector or NULL] Coordinate of lower corner of one
	#'        cell. If NULL, c(0,...,0) is used
	bin_origin = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description 
    #' Create a new OTHist object.
	#' @param p [vector] Vector of probability
	#' @param c [matrix] Vector of center of bins, with nrow = n_samples and ncol = n_features
	#'
	#' @return A new `OTHist` object.
	initialize = function( p , c )
	{
		self$p = p
		self$c = c
	}
	##}}}
	
	)
	
	
	######################
	## Private elements ##
	######################
)
