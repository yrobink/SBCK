
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
#' This base class can be considered as the identity pre-post processing, and
#' is used to be herited by others pre/post processing class. The key ideas are:\cr
#' - A PrePostProcessing based class contains a bias correction method, initalized
#'   by the `bc_method` argument, always available for all herited class\cr
#' - The `pipe` keyword is a list of pre/post processing class, applied one after
#'   the other.\cr
#' 
#' Try with an example, start with a dataset similar to tas/pr:\cr
#' >>> XY = SBCK::dataset_like_tas_pr(2000)\cr
#' >>> X0 = XY$X0\cr
#' >>> X1 = XY$X1\cr
#' >>> Y0 = XY$Y0\cr
#' 
#' The first column is Gaussian, but the second is an exponential law with a Dirac
#' mass at 0, represented the 0 of precipitations. For a quantile mapping
#' correction in the calibration period, we just apply\cr
#' >>> qm = SBCK::QM$new()\cr
#' >>> qm$fit(Y0,X0)\cr
#' >>> Z0 = qm$predict(X0)\cr
#' 
#' Now, if we want to pre-post process with the SSR method (0 are replaced by
#' random values between 0 (excluded) and the minimal non zero value), we write:\cr
#' >>> ppp = SBCK::PPPSSR$new( bc_method = QM , cols = 2 )\cr
#' >>> ppp$fit(Y0,X0)\cr
#' >>> Z0 = ppp$predict(X0)\cr
#' 
#' The SSR approach is applied only on the second column (the precipitation), and
#' the syntax is the same than for a simple bias correction method.\cr
#' 
#' Imagine now that we want to apply the SSR, and to ensure the positivity of CDFt
#' for precipitation, we also want to use the LogLinLink pre-post processing
#' method. This can be done with the following syntax:\cr
#' >>> ppp = PPPLogLinLink$new( bc_method = CDFt , cols = 2 ,\cr
#' >>>                          pipe = list(PPPSSR) , \cr
#' >>>                          pipe_kwargs = list( list(cols = 2) ) )\cr
#' >>> ppp$fit(Y0,X0,X1)\cr
#' >>> Z = ppp$predict(X1,X0)\cr
#' 
#' With this syntax, the pre processing operation is
#' PPPLogLinLink$transform(PPPSSR$transform(data)) and post processing operation
#' PPPSSR$itransform(PPPLogLinLink$itransform(bc_data)). So the formula can read
#' from right to left (as the mathematical composition). Note it is equivalent
#' to define:\cr
#' >>> ppp = PrePostProcessing$new( bc_method = CDFt,\cr
#' >>>                              pipe = list(PPPLogLinLink,PPPSSR),\cr
#' >>>                              pipe_kwargs = list( list(cols=2) , list(cols=2) ) )\cr
#' 
#'
#' @examples
#' ## Start with data
#' XY = SBCK::dataset_like_tas_pr(2000)
#' X0 = XY$X0
#' X1 = XY$X1
#' Y0 = XY$Y0
#' 
#' ## Define pre/post processing method
#' ppp = PrePostProcessing$new( bc_method = CDFt,
#'                              pipe = list(PPPLogLinLink,PPPSSR),
#'                              pipe_kwargs = list( list(cols=2) , list(cols=2) ) )
#'
#' ## Bias correction
#' ppp$fit(Y0,X0,X1)
#' Z = ppp$predict(X1,X0)
#' 
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
			Z01t = private$.bc_method$predict(X1t,X0t)
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

