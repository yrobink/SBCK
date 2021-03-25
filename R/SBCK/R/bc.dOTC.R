
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

#' dOTC (dynamical Optimal Transport Correction) method
#'
#' @description
#' Perform a multivariate (non stationary) bias correction.
#'
#' @details
#' Three random variables are needed, Y0, X0 and X1. The dynamic between X0 and
#' X1 is estimated, and applied to Y0 to estimate Y1. Finally, OTC is used
#' between X1 and the Y1 estimated.
#'
#' @references Robin, Y., Vrac, M., Naveau, P., Yiou, P.: Multivariate
#'             stochastic bias corrections with optimal transport, Hydrol. Earth
#'             Syst. Sci., 23, 773–786, 2019,
#'             https://doi.org/10.5194/hess-23-773-2019
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#' ## Bin length
#' bin_width = c(0.2,0.2)
#'
#' ## Bias correction
#' ## Step 1 : construction of the class dOTC 
#' dotc = SBCK::dOTC$new( bin_width ) 
#' ## Step 2 : Fit the bias correction model
#' dotc$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction, Z is a list containing
#' ## corrections
#' Z = dotc$predict(X1,X0) 
#' Z$Z0 ## Correction in calibration period
#' Z$Z1 ## Correction in projection period
#'
#' @export
dOTC = R6::R6Class( "dOTC" ,
	inherit = SBCK::OTC,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new dOTC object.
	#' @param bin_width [vector or NULL] A vector of lengths of the cells
	#'        discretizing R^{numbers of variables}. If NULL, it is estimating
	#'        during the fit
	#' @param bin_origin [vector or NULL] Coordinate of lower corner of one
	#'        cell. If NULL, c(0,...,0) is used
	#' @param cov_factor [string or matrix] Covariance factor to correct the
	#'        dynamic transferred between X0 and Y0. For string, available
	#'        values are "std" and "cholesky"
	#' @param ot [OTSolver] Optimal Transport solver, default is the network
	#'        simplex
	#'
	#' @return A new `dOTC` object.
	initialize = function( bin_width = NULL , bin_origin = NULL , cov_factor = "std" , ot = SBCK::OTNetworkSimplex$new() )
	{
		super$initialize( bin_width , bin_origin , ot )
		private$cov_factor = cov_factor
	},
	##}}}
	
	## Fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return NULL
	fit = function( Y0 , X0 , X1 )
	{
		## Dimension and data formating
		if( !is.matrix(Y0) ) Y0 = matrix( Y0 , nrow = length(Y0) , ncol = 1 )
		if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		if( !is.matrix(X1) ) X1 = matrix( X1 , nrow = length(X1) , ncol = 1 )
		
		## Bin width
		if( is.null(self$bin_width) )
		{
			self$bin_width = SBCK::bin_width_estimator( list(Y0,X0,X1) )
		}
		if( is.null(self$bin_origin) )
		{
			self$bin_origin = base::rep( 0. , length(self$bin_width) )
		}
		
		self$n_features = base::ncol(Y0)
		
		## Correction factor of motion
		cf = NULL
		if( self$n_features == 1 && !is.numeric(private$cov_factor) )
		{
			cf = sd(Y0) * sd(X0)^(-1)
		}
		else
		{
			if( !is.matrix(private$cov_factor) )
			{
				if( private$cov_factor == "cholesky" )
				{
					cf = base::chol( stats::cov(Y0) ) %*% base::solve( base::chol( stats::cov(X0) ) )
				}
				else
				{
					cf = base::diag( base::apply( Y0 , 2 , sd ) * base::apply( X0 , 2 , stats::sd )^(-1) )
				}
			}
		}
		## Motion
		private$otcB0R0 = SBCK::OTC$new( self$bin_width , self$bin_origin ) ## Bias
		otcB0B1 = SBCK::OTC$new( self$bin_width , self$bin_origin ) ## Evolution
		private$otcB0R0$fit( Y0 , X0 )
		otcB0B1$fit( X1 , X0 )
		
		motion = otcB0B1$predict(X0) - X0
		if( self$n_features == 1 )
		{
			motion = cf * motion
		}
		else
		{
			motion = base::t(base::apply( motion , 1 , function(x){ invisible( cf %*% x ) } ))
		}
		
		## Estimation of Y1
		Y1 = private$otcB0R0$predict( X0 ) + motion
		
		super$fit( Y1 , X1 )
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #'
	#' Note: Only the center of the bins associated to the corrected points are
	#' returned, but all corrections of the form:
	#' >> bw = dotc$bin_width / 2
	#' >> n  = base::prod(base::dim(X1))
	#' >> Z1 = dotc$predict(X1)
	#' >> Z1 = Z1 + t(matrix(stats::runif( n = n min = - bw , max = bw ) , ncol = dim(X1)[1] ))
	#' are equivalent for OTC.
	#'
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 , X0 = NULL )
	{
		Z1 = super$predict(X1)
		
		if( !is.null(X0) )
		{
			Z0 = private$otcB0R0$predict(X0)
			return( list( Z1 = Z1 , Z0 = Z0 ) )
		}
		return(Z1)
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
	
	otcB0R0 = NULL ,
	cov_factor = NULL
	)
)
