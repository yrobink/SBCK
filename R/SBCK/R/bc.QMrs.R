
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

#' Quantile Mapping RankShuffle method
#'
#' @description
#' Perform a multivariate bias correction of X with respect to Y 
#'
#' @details
#' Dependence is corrected with multi_schaake_shuffle.
#'
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
#'
#' ## Bias correction
#' ## Step 1 : construction of the class QMrs 
#' qmrs = SBCK::QMrs$new() 
#' ## Step 2 : Fit the bias correction model
#' qmrs$fit( Y0 , X0 )
#' ## Step 3 : perform the bias correction
#' Z0 = qmrs$predict(X0)
#'
#' @export
QMrs = R6::R6Class( "QMrs" ,
	
	inherit = SBCK::QM,
	
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
    #' Create a new QMrs object.
	#' @param irefs [vector of int] Indexes for shuffle. Defaults is base::c(1)
	#'        model
	#' @param ... [] all others arguments are passed to QM class.
	#'
	#' @return A new `QMrs` object.
	initialize = function( irefs = base::c(1) , ... )
	{
		super$initialize(...)
		self$irefs = irefs
		private$ssr = SBCK::SchaakeShuffleRef$new( 1 )
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #'
    #' @return NULL
	fit = function( Y0 , X0 )
	{
		super$fit( Y0 , X0 )
		private$ssr$fit(Y0)
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #'
    #' @return [matrix] Return the corrections of X0
	predict = function( X0 )
	{
		Z0u = super$predict(X0)
		Z0 = base::array( NA , dim = base::c( base::dim(Z0u) , length(self$irefs) ) )
		for( i in 1:length(self$irefs) )
		{
			private$ssr$ref = self$irefs[i]
			Z0[,,i] = private$ssr$predict(Z0u)
		}
		if( length(self$irefs) == 1 )
			return(Z0[,,1])
		return(Z0)
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
