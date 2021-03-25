
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

#' ECBC (Empirical Copula Bias Correction) method
#' 
#' @description
#' Perform a multivariate (non stationary) bias correction.
#'
#' @details
#' use Schaake shuffle
#' 
#' @references Vrac, M. and P. Friederichs, 2015: Multivariate—Intervariable,
#'             Spatial, and Temporal—Bias Correction. J. Climate, 28, 218–237,
#'             https://doi.org/10.1175/JCLI-D-14-00059.1
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
#' ## Step 1 : construction of the class ECBC
#' ecbc = SBCK::ECBC$new() 
#' ## Step 2 : Fit the bias correction model
#' ecbc$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = ecbc$predict(X1,X0) 
#'
#' @export
ECBC = R6::R6Class( "ECBC" ,
	
	inherit = SBCK::CDFt,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new ECBC object.
	#' @param ... This class is based to CDFt, and takes the same arguments.
	#' @return A new `ECBC` object.
	initialize = function(...)
	{
		kwargs = list(...)
		base::do.call( super$initialize , kwargs )
		private$ss = SBCK::SchaakeShuffle$new()
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return NULL
	fit = function( Y0 , X0 , X1 )
	{
		super$fit( Y0 , X0 , X1 )
		private$ss$fit(Y0)
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 , X0 = NULL )
	{
		Z = super$predict(X1,X0)
		if( is.null(X0) )
			return(private$ss$predict(Z))
		
		Z0 = private$ss$predict(Z$Z0)
		Z1 = private$ss$predict(Z$Z1)
		
		return( list( Z1 = Z1 , Z0 = Z0 ) )
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
	
	ss  = NULL
	)
)
