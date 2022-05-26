
## Copyright(c) 2022 Yoann Robin
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

## PPPDiffRef ##{{{

#' PPPDiffRef
#'
#' @description
#' Apply the diff w.r.t. a ref transformation.
#'
#' @details
#' Transform a dataset such that all `lower` dimensions are replaced by
#' the `ref` dimension minus the `lower`; and all `upper` dimensions are
#' replaced by `upper` minus `ref`.
#'
#' @examples
#' ## Parameters
#' size  = 2000
#' nfeat = 5
#' sign  = base::sample( base::c(-1,1) , nfeat - 1 , replace = TRUE )
#' 
#' ## Build data
#' X     = matrix( stats::rnorm( n = size ) , ncol = 1 )
#' for( s in sign )
#' {
#' 	X = base::cbind( X , X[,1] + s * base::abs(matrix( stats::rnorm(n = size) , ncol = 1 )) )
#' }
#' 
#' ## PPP
#' lower = which( sign ==  1 ) + 1
#' upper = which( sign == -1 ) + 1
#' ppp   = SBCK::PPPDiffRef$new( ref = 1 , lower = lower , upper = upper )
#' Xt    = ppp$transform(X)
#' Xti   = ppp$itransform(Xt)
#' 
#' print( base::max( base::abs( X - Xti ) ) )
#' 
#' @export
PPPDiffRef = R6::R6Class( "PPPDiffRef" ,
	
	inherit = PrePostProcessing,
	
	## Public list {{{
	
	public = list(
	
	## Public arguments
	#' @field ref [integer] The reference column
	ref   = NULL,
	#' @field lower [vector integer] Dimensions lower than ref
	lower = NULL,
	#' @field upper [vector integer] Dimensions upper than ref
	upper = NULL,
	
	## initialize ##{{{
	#' @description
    #' Create a new PPPDiffRef object.
    #' @param ref The reference column
    #' @param lower Dimensions lower than ref
    #' @param upper Dimensions upper than ref
	#' @param ... Others arguments are passed to PrePostProcessing
    #' @return A new `PPPDiffRef` object.
	initialize = function( ref , lower = NULL , upper = NULL , ... )
	{
		kwargs = list(...)
		base::do.call( super$initialize , kwargs )
		
		self$ref   = ref
		self$lower = lower
		self$upper = upper
		
		if( !is.null(self$lower) && length(self$lower) == 0 )
			self$lower = NULL
		
		if( !is.null(self$upper) && length(self$upper) == 0 )
			self$upper = NULL
	},
	##}}}
	
	## transform  ##{{{
	#' @description
    #' Apply the DiffRef transform.
    #' @param X Data to transform
    #' @return Xt a transformed matrix
	transform = function( X )
	{
		Xt = matrix( X , nrow = nrow(X) , ncol = ncol(X) )
		
		if( !is.null(self$lower) )
		{
			for( i in self$lower )
				Xt[,i] = X[,self$ref] - X[,i]
		}
		
		if( !is.null(self$upper) )
		{
			for( i in self$upper )
				Xt[,i] = X[,i] - X[,self$ref]
		}
		
		return(Xt)
	},
	##}}}
	
	## itransform  ##{{{
	#' @description
    #' Apply the DiffRef inverse transform.
    #' @param Xt Data to transform
    #' @return X a transformed matrix
	itransform = function( Xt )
	{
		X = matrix( Xt , nrow = nrow(Xt) , ncol = ncol(Xt) )
		
		if( !is.null(self$lower) )
		{
			for( i in self$lower )
				X[,i] = Xt[,self$ref] - Xt[,i]
		}
		
		if( !is.null(self$upper) )
		{
			for( i in self$upper )
				X[,i] = Xt[,i] + Xt[,self$ref]
		}
		
		return(X)
	}
	##}}}
	
	),
	##}}}
	
	## Private list ##{{{
	
	private = list(
	
	## Private arguments ##{{{
	
	##}}}
	
	## Private methods ##{{{
	
	
	##}}}
	
	)
	##}}}
	
)
##}}}

