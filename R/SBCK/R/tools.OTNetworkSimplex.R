
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


#' Optimal Transport Network Simplex solver
#'
#' @description
#' Solve the optimal transport problem 
#'
#' @details
#' use the network simplex algorithm
#'
#' @references Bazaraa, M. S., Jarvis, J. J., and Sherali, H. D.: Linear
#'             Programming and Network Flows, 4th edn., John Wiley & Sons, 2009.
#'
#' @examples
#' ## Define two dataset
#' X = stats::rnorm(2000)
#' Y = stats::rnorm(2000 , mean = 5 )
#' bw = base::c(0.1)
#' muX = SBCK::SparseHist( X , bw )
#' muY = SBCK::SparseHist( Y , bw )
#' 
#' ## Find solution
#' ot = OTNetworkSimplex$new()
#' ot$fit( muX , muY )
#' 
#' print( sum(ot$plan) ) ## Must be equal to 1
#' print( ot$success )   ## If solve is success
#' print( sqrt(sum(ot$plan * ot$C)) ) ## Cost of plan
#'
#' @export
OTNetworkSimplex = R6::R6Class( "OTNetworkSimplex" ,
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field p [double] Power of the plan
	p       = NULL,
	#' @field plan [matrix] transport plan
	plan    = NULL,
	#' @field success [bool] If the fit is a success or not
	success = NULL,
	#' @field C [matrix] Cost matrix
	C       = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description 
    #' Create a new OTNetworkSimplex object.
	#' @param p [double] Power of the plan
	#'
	#' @return A new `OTNetworkSimplex` object.
	initialize = function( p = 2 )
	{
		self$p = p
	},
	##}}}
	
	## fit ##{{{
	#' @description
	#' Fit the OT plan
	#' @param muX0 [SparseHist or OTHist] Source histogram to move
	#' @param muX1 [SparseHist or OTHist] Target histogram
	#' @param C [matrix or NULL] Cost matrix (without power p) between muX0 and
	#'       muX1, if NULL pairwise_distances is called with Euclidean distance.
	#'
	#' @return NULL
	fit = function( muX0 , muX1 , C = NULL )
	{
		self$C = C
		if( is.null(self$C) )
			self$C = SBCK::pairwise_distances( muX0$c , muX1$c )^self$p
		out = SBCK::network_simplex( muX0$p , muX1$p , self$C )
		self$plan    = out$plan
		self$success = out$success == TRUE
	}
	##}}}
	
	)
	
	
	######################
	## Private elements ##
	######################
	
)
