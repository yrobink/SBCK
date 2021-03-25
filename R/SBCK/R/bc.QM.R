
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


#' Quantile Mapping method
#'
#' @description
#' Perform an univariate bias correction of X0 with respect to Y0 
#'
#' @details
#' Correction is applied margins by margins.
#'
#' @references Panofsky, H. A. and Brier, G. W.: Some applications of statistics
#'             to meteorology, Mineral Industries Extension Services, College of
#'             Mineral Industries, Pennsylvania State University, 103 pp., 1958.
#'
#' @references Wood, A. W., Leung, L. R., Sridhar, V., and Lettenmaier, D. P.:
#'             Hydrologic Implications of Dynamical and Statistical Approaches
#'             to Downscaling Climate Model Outputs, Clim. Change, 62, 189–216,
#'             https://doi.org/10.1023/B:CLIM.0000013685.99609.9e, 2004.
#'
#' @references Déqué, M.: Frequency of precipitation and temperature extremes
#'             over France in an anthropogenic scenario: Model results and
#'             statistical correction according to observed values, Global
#'             Planet. Change, 57, 16–26,
#'             https://doi.org/10.1016/j.gloplacha.2006.11.030, 2007.
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#'
#' ## Bias correction
#' ## Step 1 : construction of the class QM 
#' qm = SBCK::QM$new() 
#' ## Step 2 : Fit the bias correction model
#' qm$fit( Y0 , X0 )
#' ## Step 3 : perform the bias correction, Z0 is the correction of
#' ## X0 with respect to the estimation of Y0
#' Z0 = qm$predict(X0)
#'
#' # ## But in fact the laws are known, we can fit parameters:
#' distY0 = list( ROOPSD::Exponential , ROOPSD::Normal )
#' distX0 = list( ROOPSD::Normal , ROOPSD::Exponential )
#' qm_fix = SBCK::QM$new( distY0 = distY0 , distX0 = distX0 )
#' qm_fix$fit( Y0 , X0 )
#' Z0 = qm_fix$predict(X0) 
#' @export
QM = R6::R6Class( "QM",
	## Public elements
	##============={{{
	
	public = list(
	
	## Arguments
	##==========
	#' @field distX0 [ROOPSD distribution or a list of them] Describe the law of
	#'        each margins. A list permit to use different laws for each
	#'        margins. Default is ROOPSD::rv_histogram.
	distX0     = NULL,
	#' @field distY0 [ROOPSD distribution or a list of them] Describe the law of
	#'        each margins. A list permit to use different laws for each
	#'        margins. Default is ROOPSD::rv_histogram.
	distY0     = NULL,
	#' @field n_features [integer] Numbers of features
	n_features = NULL,
	#' @field tol [double] Floatting point tolerance
	tol        = 1e-3,
	
	## Constructor
	##============
	
	## initialize ##{{{
	#' @description
    #' Create a new QM object.
	#' @param distX0 [ROOPSD distribution or a list of them] Describe the law of
	#'        model
	#' @param distY0 [ROOPSD distribution or a list of them] Describe the law of
	#'        observations
	#' @param ... [] kwargsX0 or kwargsY0, arguments passed to distX0 and distY0
	#'
	#' @return A new `QM` object.
	initialize = function( distX0 = ROOPSD::rv_histogram , distY0 = ROOPSD::rv_histogram , ... )
	{
		kwargs = list(...)
		self$distX0 = DistHelper$new( dist = distX0 , kwargs = kwargs[["kwargsX0"]] )
		self$distY0 = DistHelper$new( dist = distY0 , kwargs = kwargs[["kwargsY0"]] )
		self$n_features = kwargs[["n_features"]]
		self$tol = if( is.null(kwargs[["tol"]]) ) 1e-3 else kwargs[["tol"]]
	},
	##}}}
	
	
	## Methods
	##========
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #'
    #' @return NULL
	fit = function( Y0 = NULL , X0 = NULL ) 
	{
		## Data in matrix
		if( !is.null(Y0) && !is.matrix(Y0) ) Y0 = base::matrix( Y0 , ncol = 1 , nrow = length(Y0) )
		if( !is.null(X0) && !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		
		## Find n_features
		if( is.null(self$n_features ) )
		{
			if( !is.null(Y0) )
			{
				self$n_features = base::ncol(Y0)
			}
			else if( !is.null(X0) )
			{
				self$n_features = base::ncol(X0)
			}
			else
			{
				base::stop( "QM fit: if X0 = Y0 = NULL, n_features must be set during intialization" )
			}
		}
		
		## distX and distY into list
		self$distY0$set_features(self$n_features)
		self$distX0$set_features(self$n_features)
		
		## Now fit itself
		for( i in 1:self$n_features )
		{
			self$distY0$fit( Y0[,i] , i )
			self$distX0$fit( X0[,i] , i )
		}
		
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
		if( !is.null(X0) && !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		Z0 = base::matrix( NA , nrow = base::nrow(X0) , ncol = base::ncol(X0) )
		for( i in 1:self$n_features )
		{
			cdf = self$distX0$law[[i]]$cdf( X0[,i] )
			cdf[!(cdf < 1)] = 1-self$tol
			cdf[!(cdf > 0)] = self$tol
			Z0[,i] = self$distY0$law[[i]]$icdf( cdf )
		}
		return(Z0)
	}
	##}}}
	
	),
	##}}}
	
	## Private elements
	##==============={{{
	private = list(
	)
	##}}}

)

