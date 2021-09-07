
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


#' TSMBC (Time Shifted Multivariate Bias Correction)
#'
#' @description
#' Perform a bias correction of auto-correlation
#'
#' @details
#' Correct auto-correlation with a shift approach.
#'
#' @importFrom R6 R6Class
#'
#' @references Robin, Y. and Vrac, M.: Is time a variable like the others in
#'             multivariate statistical downscaling and bias correction?, Earth
#'             Syst. Dynam. Discuss. [preprint],
#'             https://doi.org/10.5194/esd-2021-12, in review, 2021.
#'
#' @examples
#' 
#' ## arima model parameters
#' modelX0 = list( ar = base::c(  0.6 , 0.2 , -0.1 ) )
#' modelY0 = list( ar = base::c( -0.3 , 0.4 , -0.2 ) )
#'
#' ## arima random generator
#' rand.genX0 = function(n){ return(stats::rnorm( n , mean = 0.2 , sd = 1   )) }
#' rand.genY0 = function(n){ return(stats::rnorm( n , mean = 0   , sd = 0.7 )) }
#'
#' ## Generate two AR processes
#' X0 = stats::arima.sim( n = 2000 , model = modelX0 , rand.gen = rand.genX0 )
#' Y0 = stats::arima.sim( n = 2000 , model = modelY0 , rand.gen = rand.genY0 )
#' X0 = as.vector( X0 )
#' Y0 = as.vector( Y0 + 5 )
#' 
#' ## And correct it with 30 lags
#' tsbc = SBCK::TSMBC$new( 30 )
#' tsbc$fit( Y0 , X0 )
#' Z0 = tsbc$predict(X0)
#'
#' @export
TSMBC = R6::R6Class( "TSMBC" ,
	
	active = list(
	
	## method ##{{{
	method = function(value)
	{
		if(missing(value))
			return(self$shift$method)
		else
			self$shift$method = value
	},
	##}}}
	
	## ref ##{{{
	ref = function(value)
	{
		if(missing(value))
			return(self$shift$ref)
		else
			self$shift$ref = value
	}
	##}}}
	
	),
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field shift [Shift class] Shift class to shift data.
	shift     = NULL,
	#' @field bc_method [SBCK::BC_method] Underlying bias correction method.
	bc_method = NULL,
	#' @field method [character] If inverse is by row or column, see class Shift
	#' @field ref [integer] reference column/row to inverse shift, see class
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new TSMBC object.
	#' @param lag [integer] max lag of autocorrelation
	#' @param bc_method [SBCK::BC_METHOD] bias correction method to use after 
	#'        shift of data, default is OTC
	#' @param method [character] If inverse is by row or column, see class Shift
	#' @param ref [integer] reference column/row to inverse shift, see class
	#'        Shift. Default is 0.5 * (lag+1)
	#' @param ... [] All others arguments are passed to bc_method
	#'
	#' @return A new `TSMBC` object.
	initialize = function( lag , bc_method = OTC , method = "row" , ref = "middle" , ... )
	{
		bc_method_args = list(...)
		self$bc_method = base::do.call( bc_method$new , bc_method_args )
		if( ref == "middle" )
			ref = as.integer(0.5 * (lag+1) )
		self$shift = Shift$new( lag , method , ref )
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
		Xs = self$shift$transform(X0)
		Ys = self$shift$transform(Y0)
		self$bc_method$fit( Ys , Xs )
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #'
    #' @return [matrix] Return the corrections of X0
	predict = function(X0)
	{
		return(self$shift$inverse( self$bc_method$predict( self$shift$transform(X0) ) )) 
	}
	##}}}
	
	)
)
