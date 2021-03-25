
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


#' OTC (Optimal Transport Correction) method
#'
#' @description
#' Perform a multivariate bias correction of X0 with respect to Y0.
#'
#' @details
#' Joint distribution, i.e. all dependence are corrected.
#'
#' @references Robin, Y., Vrac, M., Naveau, P., Yiou, P.: Multivariate
#'             stochastic bias corrections with optimal transport, Hydrol. Earth
#'             Syst. Sci., 23, 773–786, 2019,
#'             https://doi.org/10.5194/hess-23-773-2019
#'
#' @examples
#' ## Two bivariate random variables (rnorm and rexp are inverted between ref 
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#'
#' ## Bin length
#' bin_width = SBCK::bin_width_estimator( list(X0,Y0) )
#'
#' ## Bias correction
#' ## Step 1 : construction of the class OTC 
#' otc = SBCK::OTC$new( bin_width ) 
#' ## Step 2 : Fit the bias correction model
#' otc$fit( Y0 , X0 )
#' ## Step 3 : perform the bias correction, Z0 is the correction of
#' ## X0 with respect to the estimation of Y0
#' Z0 = otc$predict(X0)
#'
#' @importFrom R6 R6Class
#' @export
OTC = R6::R6Class( "OTC" ,
	
	
	public = list(
	
	
	###############
	## Arguments ##
	###############
	
	#' @field bin_width [vector or NULL] A vector of lengths of the cells
	#'        discretizing R^{numbers of variables}. If NULL, it is estimating
	#'        during the fit
	bin_width  = NULL,
	#' @field bin_origin [vector or NULL] Coordinate of lower corner of one
	#'        cell. If NULL, c(0,...,0) is used
	bin_origin = NULL,
	#' @field muX [SparseHist] Histogram of the data from the model
	muX        = NULL,
	#' @field muY [SparseHist] Histogram of the data from the observations
	muY        = NULL,
	#' @field ot [OTSolver] Optimal Transport solver, default is the network
	#'        simplex
	ot         = NULL,
	#' @field plan [matrix] The plan computed by the ot solver.
	plan       = NULL,
	#' @field n_features [integer] Numbers of features
	n_features = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new OTC object.
	#' @param bin_width [vector or NULL] A vector of lengths of the cells
	#'        discretizing R^{numbers of variables}. If NULL, it is estimating
	#'        during the fit
	#' @param bin_origin [vector or NULL] Coordinate of lower corner of one
	#'        cell. If NULL, c(0,...,0) is used
	#' @param ot [OTSolver] Optimal Transport solver, default is the network
	#'        simplex
	#'
	#' @return A new `OTC` object.
	initialize = function( bin_width = NULL , bin_origin = NULL , ot = SBCK::OTNetworkSimplex$new() )
	{
		self$bin_width  = bin_width
		self$bin_origin = bin_origin
		self$ot = ot
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
		## Dimension and data formating
		if( !is.matrix(Y0) ) Y0 = matrix( Y0 , nrow = length(Y0) , ncol = 1 )
		if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		self$n_features = base::ncol(Y0)
		
		if( is.null(self$bin_width) )
		{
			self$bin_width = SBCK::bin_width_estimator( list(X0,Y0) )
		}
		if( is.null(self$bin_origin) )
		{
			self$bin_origin = base::rep( 0. , length(self$bin_width) )
		}
		
		self$muX = SBCK::SparseHist( X0 , self$bin_width , self$bin_origin )
		self$muY = SBCK::SparseHist( Y0 , self$bin_width , self$bin_origin )
		self$ot$fit( self$muX , self$muY )
		
		self$plan = base::t(base::apply( base::t(self$ot$plan) , 2 , function(x) { x / base::sum(x) } ))
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
	#'
	#' Note: Only the center of the bins associated to the corrected points are
	#' returned, but all corrections of the form:
	#' >> bw = otc$bin_width / 2
	#' >> n  = base::prod(base::dim(X0))
	#' >> Z0 = otc$predict(X0)
	#' >> Z0 = Z0 + t(matrix(stats::runif( n = n min = - bw , max = bw ) , ncol = dim(X0)[1] ))
	#' are equivalent for OTC.
	#'
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #'
    #' @return [matrix] Return the corrections of X0
	predict = function( X0 )
	{
		if( !is.matrix(X0) )
			X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		
		arg_X = matrix( self$muX$argwhere(X0) , nrow = base::nrow(X0) , ncol = 1 )
		n_samples = self$muY$n_samples
		arg_Y = base::apply( arg_X , 1 , function(x) { base::sample( 1:n_samples , 1 , prob = self$plan[x,] ) } )
		
		return(self$muY$c[arg_Y,]) 
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
	
	)
)
