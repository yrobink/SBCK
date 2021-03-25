
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

#' R2D2 (Rank Resampling for Distributions and Dependences) method
#'
#' @description
#' Perform a multivariate (non stationary) bias correction.
#'
#' @details
#' Use rankshuffle in calibration and projection period with CDFt
#'
#' @references Vrac, M.: Multivariate bias adjustment of high-dimensional
#'             climate simulations: the Rank Resampling for Distributions and
#'             Dependences (R2 D2 ) bias correction, Hydrol. Earth Syst. Sci.,
#'             22, 3175–3196, https://doi.org/10.5194/hess-22-3175-2018, 2018.
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
#' ## Step 1 : construction of the class R2D2 
#' r2d2 = SBCK::R2D2$new() 
#' ## Step 2 : Fit the bias correction model
#' r2d2$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = r2d2$predict(X1,X0) 
#'
#' @export
R2D2 = R6::R6Class( "R2D2" ,
	
	inherit = SBCK::CDFt,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field irefs [vector of int] Indexes for shuffle. Defaults is base::c(1)
	irefs = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new R2D2 object.
	#' @param irefs [vector of int] Indexes for shuffle. Defaults is base::c(1)
	#'        model
	#' @param ... [] all others arguments are passed to CDFt class.
	#'
	#' @return A new `R2D2` object.
	initialize = function( irefs = base::c(1) , ... )
	{
		kwargs = list(...)
		base::do.call( super$initialize , kwargs )
		self$irefs = irefs
		private$ssr = SBCK::SchaakeShuffleRef$new( 1 )
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
		private$ssr$fit(Y0)
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
		{
			Z1 = base::array( NA , dim = base::c( base::dim(Z) , length(self$irefs) ) )
			for( i in 1:length(self$irefs) )
			{
				private$ssr$ref = self$irefs[i]
				Z1[,,i] = private$ssr$predict(Z)
			}
			if( length(self$irefs) == 1 )
				return(Z1[,,1])
			return(Z1)
		}
		
		Z0 = base::array( NA , dim = base::c( base::dim(Z$Z0) , length(self$irefs) ) )
		Z1 = base::array( NA , dim = base::c( base::dim(Z$Z1) , length(self$irefs) ) )
		for( i in 1:length(self$irefs) )
		{
			private$ssr$ref = self$irefs[i]
			Z0[,,i] = private$ssr$predict(Z$Z0)
			Z1[,,i] = private$ssr$predict(Z$Z1)
		}
		if( length(self$irefs) == 1 )
		{
			Z0 = Z0[,,1]
			Z1 = Z1[,,1]
		}
		
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
	
	ssr  = NULL
	)
)
