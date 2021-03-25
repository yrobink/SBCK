
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

#' RBC (Random Bias Correction) method
#'
#' @description
#' Perform a multivariate bias correction of X with respect to Y randomly.
#'
#' @details
#' Only for comparison.
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
#' rbc = SBCK::RBC$new() 
#' ## Step 2 : Fit the bias correction model
#' rbc$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = rbc$predict(X1,X0) 
#' ## Z$Z0 # BC of X0
#' ## Z$Z1 # BC of X1
#' @export
RBC = R6::R6Class( "RBC" ,
	
	
	public = list(
	
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new RBC object.
	#'
	#' @return A new `RBC` object.
	initialize = function()
	{},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection, can be
    #'        NULL for stationary BC method
    #' @return NULL
	fit = function( Y0 , X0 , X1 = NULL )
	{
		private$Y0 = if( !is.matrix(Y0) ) matrix( Y0 , nrow = length(Y0) , ncol = 1 ) else Y0
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction. Use named keywords to use stationary or
    #' non-stationary method.
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features or NULL] Model in projection
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 = NULL , X0 = NULL )
	{
		Z0 = NULL
		Z1 = NULL
		if( !is.null(X0) )
		{
			if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
			idx = base::sample( 1:base::nrow(private$Y0) , base::nrow(X0) , replace = TRUE )
			Z0 = private$Y0[idx,]
		}
		
		if( !is.null(X1) )
		{
			if( !is.matrix(X1) ) X1 = matrix( X1 , nrow = length(X0) , ncol = 1 )
			idx = base::sample( 1:base::nrow(private$Y0) , base::nrow(X1) , replace = TRUE )
			Z1 = private$Y0[idx,]
		}
		
		if( !is.null(X0) && !is.null(X1) )
		{
			return( list( Z0 = Z0 , Z1 = Z1 ) )
		}
		else if( is.null(X0) )
		{
			return(Z1)
		}
		else
		{
			return(Z0)
		}
	}
	##}}}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	Y0 = NULL
	
	)
)
