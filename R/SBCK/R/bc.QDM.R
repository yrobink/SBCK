
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

#' QDM (Quantile delta mapping method)
#'
#' @description
#' Perform a bias correction.
#'
#' @details
#' Mix of delta and quantile method
#'
#' @references Cannon, A. J., Sobie, S. R., and Murdock, T. Q.: Bias correction
#'             of simulated precipitation by quantile mapping: how well do
#'             methods preserve relative changes in quantiles and extremes?, J.
#'             Climate, 28, 6938–6959,
#'             https://doi.org/10.1175/JCLI-D-14- 00754.1, 2015.
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#' ## Bias correction
#' ## Step 1 : construction of the class QDM
#' qdm = SBCK::QDM$new() 
#' ## Step 2 : Fit the bias correction model
#' qdm$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction, Z is a list containing
#' ## corrections
#' Z = qdm$predict(X1,X0) 
#' Z$Z0 ## Correction in calibration period
#' Z$Z1 ## Correction in projection period
#'
#' @export
QDM = R6::R6Class( "QDM" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new QDM object.
	#' @param delta [character or list] If character : "additive" or
	#'        "multiplicative". If a list is given, delta[[1]] is the delta
	#'        transform operator, and delta[[2]] its inverse.
	#' @param ... [] Named arguments passed to quantile mapping
	#'
	#' @return A new `QDM` object.
	initialize = function( delta = "additive" , ... )
	{
		## Initialize delta method
		if( class(delta) == "list" )
		{
			private$delta_method = delta[[1]]
			private$idelta_method = delta[[2]]
		}
		else if( delta == "multiplicative" )
		{
			private$delta_method  = private$mult
			private$idelta_method = private$div
		}
		else
		{
			private$delta_method  = private$add
			private$idelta_method = private$sub
		}
		
		private$qm_args = list(...)
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #'
    #' @return NULL
	fit = function( Y0 , X0 , X1 )
	{
		if( !is.matrix(Y0) ) Y0 = base::matrix( Y0 , ncol = 1 , nrow = length(Y0) )
		if( !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		if( !is.matrix(X1) ) X1 = base::matrix( X1 , ncol = 1 , nrow = length(X1) )
		
		## Fit calibration part
		private$qmX0Y0 = base::do.call( QM$new , private$qm_args )
		private$qmX0Y0$fit( Y0 , X0 )
		
		## Fit delta
		qmX1X0 = base::do.call( QM$new , private$qm_args )
		qmX1X0$fit( X0 , X1 )
		private$delta = private$idelta_method( X1 , qmX1X0$predict(X1) )
		
		## Fit projection part
		private$qmX1Y0 = base::do.call( QM$new , private$qm_args )
		private$qmX1Y0$fit( Y0 , X1 )
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #'
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 , X0 = NULL )
	{
		if( !is.null(X0) && !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		if( !is.matrix(X1) ) X1 = base::matrix( X1 , ncol = 1 , nrow = length(X1) )
		
		Z1 = private$delta_method( private$qmX1Y0$predict(X1) , private$delta )
		if( !is.null(X0) )
		{
			Z0 = private$qmX0Y0$predict(X0)
			return( list( Z1 = Z1 , Z0 = Z0 ) )
		}
		return(Z1)
	}
	##}}}
	
	),
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	delta_method  = NULL,
	idelta_method = NULL,
	delta         = NULL,
	qm_args       = NULL,
	qmX0Y0        = NULL,
	qmX1Y0        = NULL,
	
	
	#############
	## Methods ##
	#############
	
	add = function(x,y)##{{{
	{ 
		return( x + y )
	},
	##}}}
	
	mult = function(x,y)##{{{
	{ 
		return( x * y )
	},
	##}}}
	
	sub = function(x,y)##{{{
	{ 
		return( x - y )
	},
	##}}}
	
	div = function(x,y)##{{{
	{ 
		return( x / y )
	}
	##}}}
	
	)
)
