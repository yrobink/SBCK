
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

#' AR2D2 (Analogues Rank Resampling for Distributions and Dependences) method
#' 
#' @description
#' Perform a multivariate (non stationary) bias correction.
#'
#' @details
#' Use Quantiles shuffle in calibration and projection period with CDFt
#'
#' @references Vrac, M. et S. Thao (2020). “R2 D2 v2.0 : accounting for temporal
#'             dependences in multivariate bias correction via analogue rank
#'             resampling”. In : Geosci. Model Dev. 13.11, p. 5367-5387.
#'             doi :10.5194/gmd-13-5367-2020.
#'
#' @examples
#' 
#' ## Three 4-variate random variables
#' Y0 = matrix( stats::rnorm( n = 1000 ) , ncol = 4 ) ## Biased in calibration period
#' X0 = matrix( stats::rnorm( n = 1000 ) , ncol = 4 ) / 2 + 3 ## Reference in calibration period
#' X1 = matrix( stats::rnorm( n = 1000 ) , ncol = 4 ) * 2 + 6 ## Biased in projection period
#'
#' ## Bias correction
#' cond_col = base::c(2,4)
#' lag_search = 6
#' lag_keep = 3
#' ## Step 1 : construction of the class AR2D2 
#' ar2d2 = SBCK::AR2D2$new( cond_col , lag_search , lag_keep ) 
#' ## Step 2 : Fit the bias correction model
#' ar2d2$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = ar2d2$predict(X1,X0) 
#'
#' @export
AR2D2 = R6::R6Class( "AR2D2" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field mvq [MVQuantilesShuffle] Class to transform dependance structure
	mvq       = NULL,
	#' @field bc_method [SBCK::] Bias correction method
	bc_method = NULL,
	#' @field bckwargs [list] List of arguments of bias correction
	bckwargs  = NULL,
	#' @field bcm_ [SBCK::] Instancied bias correction method
	bcm_      = NULL,
	#' @field reverse [bool] If we apply bc_method first and then shuffle, or reverse
	reverse   = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new AR2D2 object.
    #' @param col_cond Conditionning colum
	#' @param lag_search  Number of lags to transform the dependence structure
	#' @param lag_keep Number of lags to keep
	#' @param bc_method Bias correction method
	#' @param shuffle Shuffle method used, can be quantile or rank
	#' @param reverse If we apply bc_method first and then shuffle, or reverse
	#' @param ... Others named arguments passed to bc_method$new
    #' @return A new `AR2D2` object.
	initialize = function( col_cond = base::c(1) , lag_search = 1 , lag_keep = 1 , bc_method = SBCK::CDFt , shuffle = "quantile" , reverse = FALSE , ... ) 
	{
		if( shuffle == "quantile" )
			self$mvq = MVQuantilesShuffle$new( col_cond , lag_search , lag_keep )
		else
			self$mvq = MVRanksShuffle$new( col_cond , lag_search , lag_keep )
		self$bc_method = bc_method
		self$bckwargs = list(...)
		self$reverse = reverse
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method. If X1 is NULL, the method is considered
    #' as stationary
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return NULL
	fit = function( Y0 , X0 , X1 = NULL )
	{
		self$mvq$fit(Y0)
		self$bcm_ = do.call( self$bc_method$new , self$bckwargs )
		if( is.null(X1) )
		{
			if( self$reverse )
			{
				Z0 = self$mvq$transform(X0)
				self$bcm_$fit( Y0 , Z0 )
			}
			else
			{
				self$bcm_$fit( Y0 , X0 )
			}
		}
		else
		{
			if( self$reverse )
			{
				Z0 = self$mvq$transform(X0)
				Z1 = self$mvq$transform(X1)
				self$bcm_$fit( Y0 , Z0 , Z1 )
			}
			else
			{
				self$bcm_$fit( Y0 , X0 , X1 )
			}
		}
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features or NULL] Model in projection
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL (and vice-versa), else return a list
    #'                          containing Z1 and Z0, the corrections of X1 and X0
	predict = function( X1 = NULL , X0 = NULL )
	{
		if( is.null(X0) && is.null(X1) )
		{
			return(NULL)
		}
		if( is.null(X0) )
		{
			if( self$reverse )
			{
				Z1 = self$mvq$transform(X1)
				return( self$bcm_$predict(Z1) )
			}
			else
			{
				Z1b = self$bcm_$predict(X1)
				return( self$mvq$transform(Z1b) )
			}
		}
		if( is.null(X1) )
		{
			if( self$reverse )
			{
				Z0 = self$mvq$transform(X0)
				return( self$bcm_$predict(Z0) )
			}
			else
			{
				Z0b = self$bcm_$predict(X0)
				return( self$mvq$transform(Z0b) )
			}
		}
		
		if( self$reverse )
		{
			Z0 = self$mvq$transform(X0)
			Z1 = self$mvq$transform(X1)
			return( self$bcm_$predict(Z1,Z0) )
		}
		else
		{
			Zb = self$bcm_$predict(X1,X0)
			Z  = list( Z1 = self$mvq$transform(Zb$Z1) , Z0 = self$mvq$transform(Zb$Z0) )
			return(Z)
		}
	}
	##}}}
	
	)
	
)

