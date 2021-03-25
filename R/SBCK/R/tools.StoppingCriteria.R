
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


#' Slope stopping criteria
#'
#' @description
#' Class which send a stop signal when a time series stay constant.
#'
#' @details
#' Test the slope.
#'
#' @examples
#' stop_slope = SlopeStoppingCriteria$new( 20 , 500 , 1e-3 )
#' x = 0
#' while(!stop_slope$stop)
#' {
#' 	stop_slope$append(base::exp(-x))
#' 	x = x + 0.1
#' }
#' print(stop_slope$nit)
#' @export
SlopeStoppingCriteria = R6::R6Class( "SlopeStoppingCriteria" ,
	
	public = list(
	
	################
	## Parameters ##
	################
	
	#' @field minit [integer] Minimal number of iterations. At least 3.
	minit    = NULL,
	#' @field maxit [integer] Maximal number of iterations.
	maxit    = NULL,
	#' @field nit [integer] Number of iterations.
	nit      = NULL,
	#' @field tol [float] Tolerance to control if slope is close to zero
	tol      = NULL,
	#' @field stop [bool] If we stop
	stop     = NULL,
	#' @field criteria [vector] State of criteria
	criteria = NULL,
	#' @field slope [vector] Values of slope
	slope    = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new SlopeStoppingCriteria object.
	#' @param minit [integer] Minimal number of iterations. At least 3.
	#' @param maxit [integer] Maximal number of iterations.
	#' @param tol [float] Tolerance to control if slope is close to zero
	#'
	#' @return A new `SlopeStoppingCriteria` object.
	initialize = function( minit , maxit , tol )
	{
		self$minit = max(minit,3)
		self$maxit = maxit
		self$tol = tol
		self$nit = 1
		self$stop = FALSE
		self$criteria = numeric(maxit)
		self$slope = numeric(maxit)
	},
	##}}}
	
	## reset ##{{{
	#' @description
    #' Reset the class
	#'
	#' @return NULL
	reset = function()
	{
		self$nit = 1
		self$stop = FALSE
		self$criteria = numeric(self$maxit)
		self$slope = numeric(self$maxit)
	},
	##}}}
	
	## append ##{{{
	#' @description
    #' Add a new value
    #' @param value [double] New metrics
	#'
	#' @return NULL
	append = function(value)
	{
		self$criteria[self$nit] = value
		if( self$nit >= 2 )
		{
			x = 1:self$nit
			y = self$criteria[1:self$nit]
			lm = stats::lm( y ~ x )
			self$slope[self$nit-1] = as.vector(lm$coefficients[2])
		}
		self$stop = (self$nit >= self$minit) && (self$nit >= self$maxit - 1 || base::abs(self$slope[self$nit-1]) < self$tol)
		self$nit = self$nit + 1
	}
	##}}}
	
	)
)


