
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

## SchaakeShuffle ##{{{


#' ShaakeShuffle class
#'
#' @description
#' Perform the Schaake Shuffle 
#'
#' @details
#' as fit/predict mode
#'
#'
#' @examples
#' X0 = matrix( stats::runif(20) , ncol = 2 )
#' Y0 = matrix( stats::runif(20) , ncol = 2 )
#' ss = SchaakeShuffle$new()
#' ss$fit(Y0)
#' Z0 = ss$predict(X0)
#'
#' @export
SchaakeShuffle = R6::R6Class( "SchaakeShuffle" ,
	
	public = list( ##{{{
	
	## initialize ##{{{
	#' @description
    #' Create a new ShaakeShuffle object.
	#' @param Y0 [vector] The reference vector
	#'
	#' @return A new `ShaaleShuffle` object.
	initialize = function( Y0 = NULL )
	{
		private$Y0 = NULL
		if( !is.null(Y0) )
			self$fit(Y0)
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the model
	#' @param Y0 [vector] The reference vector
	#'
	#' @return NULL
	fit = function( Y0 )
	{
		private$Y0 = Y0
		if( !is.matrix(Y0) ) private$Y0 = matrix( Y0 , nrow = length(Y0) , ncol = 1 )
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Fit the model
	#' @param X0 [vector] The vector to apply shuffle
	#'
	#' @return Z0 [vector] data shuffled
	predict = function( X0 )
	{
		if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		
		YY = NULL
		XX = NULL
		nrowY0 = base::nrow(private$Y0)
		nrowX0 = base::nrow(X0)
		n_features = base::ncol(private$Y0)
		if( nrowY0 < nrowX0 )
		{
			YY = matrix( NA , nrow = nrowX0 , ncol = n_features )
			YY[1:nrowY0,] = private$Y0
			YY[nrowY0:nrowX0,] = private$Y0[base::sample( nrowY0 , nrowX0 - nrowY0 , replace = TRUE ),]
			XX = X0
		}
		else if( nrowX0 < nrowY0 )
		{
			XX = matrix( NA , nrow = nrowY0 , ncol = n_features )
			XX[1:nrowX0,] = X0
			XX[nrowX0:nrowY0,] = X0[base::sample( nrowX0 , nrowY0 - nrowX0 , replace = TRUE ),]
			YY = private$Y0
		}
		else
		{
			XX = X0
			YY = private$Y0
		}
		
		ZZ = matrix( NA , nrow = base::nrow(XX) , ncol = base::ncol(XX) )
		for( i in 1:n_features )
		{
			ZZ[,i] = private$predict_univariate( YY[,i] , XX[,i] )
		}
		Z = ZZ[1:nrowX0,]
		
		return(Z)
	}
	##}}}
	
	),
	##}}}
	
	private = list(##{{{
	
	###############
	## Arguments ##
	###############
	
	Y0 = NULL,
	
	
	#############
	## Methods ##
	#############
	
	predict_univariate = function( Y , X )##{{{
	{
		rank_X  = base::rank(X)
		rank_Y  = base::rank(Y)
		arank_X = base::order(rank_X)
		Z = X[arank_X][rank_Y]
		return(Z)
	}
	##}}}
	
	)
	##}}}
)
##}}}

## schaake_shuffle ##{{{

#' schaake_shuffle function
#'
#' Apply the Schaake shuffle to transform the rank of X0 such that its
#' correspond to the rank of Y0
#' 
#' @usage schaake_shuffle(Y0,X0)
#' @param X0 [vector] The vector to transform the rank
#' @param Y0 [vector] The reference vector
#'
#' @return Z0 [vector] X shuffled.
#'
#' @examples
#' X0 = stats::runif(10)
#' Y0 = stats::runif(10)
#' Z0 = SBCK::schaake_shuffle( Y0 , X0 )
#'
#' @export
schaake_shuffle = function( Y0 , X0 )
{
	ss = SchaakeShuffle$new(Y0)
	return(ss$predict(X0))
}
##}}}

## SchaakeShuffleRef ##{{{

#' ShaakeShuffleRef class
#'
#' @description
#' Match the rank structure of X with them of Y by reordering X.
#'
#' @details
#' Fix one features to keep the structure of X.
#'
#' @param ref [integer] The reference
#' @param X0 [vector] The vector to transform the rank
#' @param Y0 [vector] The reference vector
#'
#' @examples
#' X0 = matrix( stats::runif(20) , ncol = 2 )
#' Y0 = matrix( stats::runif(20) , ncol = 2 )
#' ss = SchaakeShuffleRef$new( ref = 1 )
#' ss$fit(Y0)
#' Z0 = ss$predict(X0)
#'
#' @export
SchaakeShuffleRef = R6::R6Class( "SchaakeShuffleRef" ,
	
	inherit = SchaakeShuffle,
	
	public = list( ##{{{
	
	#' @field ref [integer] Reference
	ref = NULL,
	
	## initialize ##{{{
	#' @description
    #' Create a new ShaakeShuffleRef object.
	#' @param ref [integer] Reference
	#' @param Y0 [vector] The reference vector
	#'
	#' @return A new `ShaaleShuffleRef` object.
	initialize = function( ref , Y0 = NULL )
	{
		self$ref = ref
		super$initialize(Y0)
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the model
	#' @param Y0 [vector] The reference vector
	#'
	#' @return NULL
	fit = function(Y0)
	{
		super$fit(Y0)
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Fit the model
	#' @param X0 [vector] The vector to apply shuffle
	#'
	#' @return Z0 [vector] data shuffled
	predict = function(X0)
	{
		if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		Z0 = super$predict(X0)
		
		rank_X0  = base::rank(X0[,self$ref])
		rank_Z0  = base::rank(Z0[,self$ref])
		arank_Z0 = base::order(rank_Z0)
		Z0 = Z0[arank_Z0,][rank_X0,]
		return(Z0)
	}
	##}}}
	
	)
	##}}}

)
##}}}

