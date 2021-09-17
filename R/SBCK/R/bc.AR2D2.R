
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

#' CDFt method (Cumulative Distribution Function transfer)
#'
#' @description
#' Perform an univariate bias correction of X with respect to Y.
#'
#' @details
#' Correction is applied margins by margins.
#'
#' @references Michelangeli, P.-A., Vrac, M., and Loukos, H.: Probabilistic 
#'             downscaling approaches: Application to wind cumulative 
#'             distribution functions, Geophys. Res. Lett., 36, L11708, 
#'             https://doi.org/10.1029/2009GL038401, 2009.
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#' ## Bias correction
#' ## Step 1 : construction of the class AR2D2
#' lag_search = 3
#' lag_keep   = 2
#' ar2d2 = SBCK::AR2D2$new(lag_search,lag_keep) 
#' ## Step 2 : Fit the bias correction model
#' ar2d2$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction, Z is a list containing
#' ## corrections
#' Z = ar2d2$predict(X1,X0) 
#' Z$Z0 ## Correction in calibration period
#' Z$Z1 ## Correction in projection period
#' @export
AR2D2 = R6::R6Class( "AR2D2" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field bc_univ [SBCK::] Bias correction method
	bc_univ = NULL,
	#' @field ssmr [SchakeShuffleMultiRef] Multivariate shuffle
	ssmr    = NULL,
	#' @field debug [list] Debug list
	debug   = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize  ##{{{
	#' @description
    #' Create a new AR2D2 object.
	#' @param lag_search An integer corresponding to the number of time lags to
	#'            account for when searching in \code{refdata} the best analogue of the
	#'            conditioning dimension rank association observed in \code{bc1d}. The default
	#'            value is no lag, i.e., \code{lag_search}=0.
	#' @param lag_keep An integer corresponding to the number of time lags to keep
	#'            in the correction for each best analogue of rank associations found.
	#'            \code{lag_keep} has to be lower or equal to \code{lag_search}.  The default
	#'            value is no lag, i.e., \code{lag_search}=0.
	#' @param cond_cols [vector of integer] Columns for conditioning
	#' @param bc_univ [BC method] A bias correction method
	#' @param ... Others arguments are passed to bc_univ
	#' @return A new `AR2D2` object.
	initialize = function( lag_search , lag_keep , cond_cols = base::c(1) , bc_univ = SBCK::CDFt , ... )
	{
		self$bc_univ = base::do.call( bc_univ$new , list(...) )
		self$ssmr    = SchaakeShuffleMultiRef$new( lag_search , lag_keep , cond_cols )
		self$debug   = list()
	},
	##}}}
	
	## fit  ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return NULL
	fit = function( Y0 , X0 , X1 )
	{
		self$bc_univ$fit(Y0,X0,X1)
		self$ssmr$fit(Y0)
	},
	##}}}
	
	## predict  ##{{{
	#' @description
    #' Predict the correction
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 , X0 = NULL )
	{
		if( is.null(X0) )
		{
			Z1u = self$bc_univ$predict(X1)
			Z1  = self$ssmr$predict(Z1u)
			self$debug$Z1u = Z1u
			return(Z1)
		}
		else
		{
			Zu = self$bc_univ$predict(X1,X0)
			Z0  = self$ssmr$predict(Zu$Z0)
			Z1  = self$ssmr$predict(Zu$Z1)
			return( list( Z1 = Z1 , Z0 = Z0 ) )
		}
	}
	##}}}
	
	)
)

