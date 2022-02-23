
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


## PrePostProcessing ##{{{

#' PrePostProcessing base class
#'
#' @description
#' Base class to pre/post process data before/after a bias correction
#'
#' @details
#' Equivalent to the identity
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#'
#' ## Bias correction
#' ## Step 1 : construction of the class RBC
#' ppp = PrePostProcessing$new( bc_method = SBCK::CDFt ) 
#' ## Step 2 : Fit the bias correction model
#' ppp$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = ppp$predict(X1,X0) 
#' ## Z$Z0 # BC of X0
#' ## Z$Z1 # BC of X1
#' @export
PrePostProcessing = R6::R6Class( "PrePostProcessing" ,
	
	## Public list {{{
	
	public = list(
	
	## initialize ##{{{
	#' @description
    #' Create a new PrePostProcessing object.
    #' @param bc_method The bias correction method
	#' @param bc_method_kwargs Dict of keyword arguments passed to bc_method
	#' @param pipe list of others PrePostProcessing class to pipe
	#' @param pipe_kwargs list of list of keyword arguments passed to each
	#'        elements of pipe
    #' @return A new `PrePostProcessing` object.
	initialize = function( bc_method = NULL , bc_method_kwargs = list() , pipe = list() , pipe_kwargs = list() )
	{
		private$.pipe = list()
		if( !is.null(bc_method) )
		{
			private$.bc_method  = base::do.call( bc_method$new , bc_method_kwargs )
		}
		if( length(pipe) > 0 )
		{
			for( i in 1:length(pipe) )
			{
				private$.pipe[[i]] = base::do.call( pipe[[i]]$new , pipe_kwargs[[i]] )
			}
		}
	},
	##}}}
	
	## transform ##{{{
	#' @description
	#' Transformation applied to data before the bias correction. Just the
	#' identity for this class
	#' @param X [matrix: n_samples * n_features]
	#' @return Xt [matrix: n_samples * n_features]
	transform = function( X )
	{
		return(X)
	},
	##}}}
	
	## itransform ##{{{
	#' @description
	#' Transformation applied to data after the bias correction. Just the
	#' identity for this class
	#' @param Xt [matrix: n_samples * n_features]
	#' @return X [matrix: n_samples * n_features]
	itransform = function( Xt )
	{
		return(Xt)
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Apply the pre processing and fit the bias correction method. If X1 is
    #' NULL, the method is considered as stationary
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return NULL
	fit = function( Y0 , X0 , X1 = NULL )
	{
		Y0t = private$pipe_transform(Y0)
		X0t = private$pipe_transform(X0)
		X1t = private$pipe_transform(X1)
		
		if( is.null(X1) )
		{
			private$.bc_method$fit( Y0t , X0t )
		}
		else
		{
			private$.bc_method$fit( Y0t , X0t , X1t )
		}
	
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction, apply pre-processing before, and post-processing
    #' after
    #' @param X1 [matrix: n_samples * n_features or NULL] Model in projection
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL (and vice-versa), else return a list
    #'                          containing Z1 and Z0, the corrections of X1 and X0
	predict = function( X1 = NULL , X0 = NULL )
	{
		X0t = private$pipe_transform(X0)
		X1t = private$pipe_transform(X1)
		Z0t = NULL
		Z1t = NULL
		
		if( is.null(X0) )
		{
			Z1t = private$.bc_method$predict(X1t)
			return( private$pipe_itransform(Z1t) )
		}
		else if( is.null(X1) )
		{
			Z0t = private$.bc_method$predict(X0t)
			return( private$pipe_itransform(Z0t) )
		}
		else
		{
			Z01t = private$.bc_method$predict(X1t,Z0t)
			Z1 = private$pipe_itransform(Z01t$Z1)
			Z0 = private$pipe_itransform(Z01t$Z0)
			return( list( Z1 = Z1 , Z0 = Z0 ) )
		}
	
	}
	##}}}
	
	),
	##}}}
	
	## Private list ##{{{
	
	private = list(
	
	## Private arguments ##{{{
	
	.pipe      = NULL,
	.bc_method = NULL,
	
	##}}}
	
	## Private methods ##{{{
	
	pipe_transform = function(X) ##{{{
	{
		if( is.null(X) )
			return(NULL)
		
		if( length(private$.pipe) == 0 )
		{
			Xt = self$transform(X)
			return(Xt)
		}
		
		Xt = X
		for( i in length(private$.pipe):1 )
			Xt = private$.pipe[[i]]$transform(Xt)
		Xt = self$transform(Xt)
		
		return(Xt)
	},
	##}}}
	
	pipe_itransform = function(Xt) ##{{{
	{
		if( is.null(Xt) )
			return(NULL)
		
		if( length(private$.pipe) == 0 )
		{
			X = self$itransform(Xt)
			return(X)
		}
		
		X = Xt
		X = self$itransform(X)
		for( i in 1:length(private$.pipe) )
			X = private$.pipe[[i]]$itransform(X)
		
		return(X)
	}
	##}}}
	
	##}}}
	
	)
	##}}}
	
)
##}}}

