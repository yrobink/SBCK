
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

################################################################################
################################################################################
##                                                                            ##
## Notes   : CDFt is the re-implementation of the function CDFt of R package  ##
##           "CDFt" developped by Mathieu Vrac, available at                  ##
##           https://cran.r-project.org/web/packages/CDFt/index.html          ##
##           This code is governed by the GNU-GPL3 license with the           ##
##           authorization of Mathieu Vrac                                    ##
##                                                                            ##
################################################################################
################################################################################

###############
## Libraries ##
###############

###############
## Functions ##
###############



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
#' ## Step 1 : construction of the class CDFt 
#' cdft = SBCK::CDFt$new() 
#' ## Step 2 : Fit the bias correction model
#' cdft$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction, Z is a list containing
#' ## corrections
#' Z = cdft$predict(X1,X0) 
#' Z$Z0 ## Correction in calibration period
#' Z$Z1 ## Correction in projection period
#' @export
CDFt = R6::R6Class( "CDFt" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field n_features [integer] Number of features
	n_features = 0,
	#' @field tol [double] Floatting point tolerance
	tol        = 1e-3,
	
	#' @field distY0 [ROOPSD distribution or a list of them] Describe the law of
	#'        each margins. A list permit to use different laws for each
	#'        margins. Default is ROOPSD::rv_histogram.
	distY0 = NULL,
	#' @field distY1 [ROOPSD distribution or a list of them] Describe the law of
	#'        each margins. A list permit to use different laws for each
	#'        margins. Default is ROOPSD::rv_histogram.
	distY1 = NULL,
	#' @field distX0 [ROOPSD distribution or a list of them] Describe the law of
	#'        each margins. A list permit to use different laws for each
	#'        margins. Default is ROOPSD::rv_histogram.
	distX0 = NULL,
	#' @field distX1 [ROOPSD distribution or a list of them] Describe the law of
	#'        each margins. A list permit to use different laws for each
	#'        margins. Default is ROOPSD::rv_histogram.
	distX1 = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new CDFt object.
	#' @param ... Optional arguments are distX0, distX1 (models in calibration
	#'            and projection period), distY0, distY1 (observations
	#'            in calibration and projection period), and kwargsX0,
	#'            ...,kwargsY1 the arguments of each respective
	#'            distribution. The type of dist* are ROOPSD distribution,
	#'            whereas kwargs* are list.
	#' @return A new `CDFt` object.
	initialize = function(...) 
	{
		kwargs = list(...)
		self$distY0 = DistHelper$new( dist = kwargs[["distY0"]] , kwargs = kwargs[["kwargsY0"]] )
		self$distY1 = DistHelper$new( dist = kwargs[["distY1"]] , kwargs = kwargs[["kwargsY1"]] )
		self$distX0 = DistHelper$new( dist = kwargs[["distX0"]] , kwargs = kwargs[["kwargsX0"]] )
		self$distX1 = DistHelper$new( dist = kwargs[["distX1"]] , kwargs = kwargs[["kwargsX1"]] )
		self$n_features = kwargs[["n_features"]]
		self$tol        = if( !is.null(kwargs[["tol"]]) ) kwargs[["tol"]] else 1e-3
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
		## Dimension and data formating
		##=============================
		if( is.vector(Y0) ) Y0 = matrix( Y0 , nrow = length(Y0) , ncol = 1 )
		if( is.vector(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		if( is.vector(X1) ) X1 = matrix( X1 , nrow = length(X1) , ncol = 1 )
		
		## Set features
		##=============
		if( is.null(self$n_features ) )
		{
			if( !is.null(Y0) )
				self$n_features = base::ncol(Y0)
			else if( !is.null(X0) )
				self$n_features = base::ncol(X0)
			else if( !is.null(X1) )
				self$n_features = base::ncol(X1)
			else
				base::stop("SBCK::CDFt: If Y0 == X0 == X1 == NULL, n_features must be set during initialization")
		}
		self$distY0$set_features(self$n_features)
		self$distY1$set_features(self$n_features)
		self$distX0$set_features(self$n_features)
		self$distX1$set_features(self$n_features)
		
		
		## Start fit itself
		##=================
		for( i in 1:self$n_features )
		{
			self$distY0$fit( Y0[,i] , i )
			self$distX0$fit( X0[,i] , i )
			self$distX1$fit( X1[,i] , i )
			
			if( self$distY1$is_frozen(i) )
			{
				self$distY1$law[[i]] = self$distY1$dist[[i]]
			}
			else
			{
				Y1 = NULL
				if( self$distY0$is_parametric(i) && self$distX0$is_parametric(i) && self$distX1$is_parametric(i) )
				{
					Y1 = self$distX1$law[[i]]$icdf( self$distX0$law[[i]]$cdf( self$distY0$law[[i]]$icdf( self$distX1$law[[i]]$cdf( X1[,i] ) ) ) )
				}
				else
				{
					Y0uni = if( is.null(Y0) ) self$distY0$law[[i]]$rvs(10000) else Y0[,i]
					X0uni = if( is.null(X0) ) self$distX0$law[[i]]$rvs(10000) else X0[,i]
					X1uni = if( is.null(X1) ) self$distX1$law[[i]]$rvs(10000) else X1[,i]
					Y1 = private$infer_Y1_new( Y0uni , X0uni , X1uni , i ) ## infer Y1 here
				}
				self$distY1$fit( Y1 , i )
			}
		}
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
		if( !is.matrix(X1) ) X1 = matrix( X1 , nrow = length(X1) , ncol = 1 )
		
		if( is.null(X0) )
		{
			Z1 = matrix( NA , nrow = base::nrow(X1) , ncol = base::ncol(X1) )
			for( i in 1:self$n_features )
			{
				cdf = self$distX1$law[[i]]$cdf( X1[,i] )
				cdf[!(cdf < 1)] = 1-self$tol
				cdf[!(cdf > 0)] = self$tol
				Z1[,i] = self$distY1$law[[i]]$icdf( cdf )
			}
			return(Z1)
		}
		else
		{
			if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
			Z0 = matrix( NA , nrow = base::nrow(X0) , ncol = base::ncol(X0) )
			Z1 = matrix( NA , nrow = base::nrow(X1) , ncol = base::ncol(X1) )
			for( i in 1:self$n_features )
			{
				cdf = self$distX0$law[[i]]$cdf( X0[,i] )
				cdf[!(cdf < 1)] = 1-self$tol
				cdf[!(cdf > 0)] = self$tol
				Z0[,i] = self$distY0$law[[i]]$icdf( cdf )
				
				cdf = self$distX1$law[[i]]$cdf( X1[,i] )
				cdf[!(cdf < 1)] = 1-self$tol
				cdf[!(cdf > 0)] = self$tol
				Z1[,i] = self$distY1$law[[i]]$icdf( cdf )
			}
			return( list( Z1 = Z1 , Z0 = Z0 ) )
		}
		
		return(Z1)
	}##}}}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(##{{{
	
	###############
	## Arguments ##
	###############
	
	qmX1Y1 = list(),
	qmX0Y0 = NULL,
	diff = NULL,
	
	
	#############
	## Methods ##
	#############
	
	infer_Y1 = function( Y0 , X0 , X1 , idx )##{{{
	{
		mY0 = base::mean(Y0)
		mX0 = base::mean(X0)
		
		X0s = X0 + mY0 - mX0
		X1s = X1 + mY0 - mX0
		
		rvY0 = self$distY0$law[[idx]]
		rvX0 = base::do.call( self$distX0$dist[[idx]]$new , self$distX0$kwargs )
		rvX0$fit(X0s)
		rvX1 = base::do.call( self$distX1$dist[[idx]]$new , self$distX1$kwargs )
		rvX1$fit(X1s)
		
		xdiff = base::abs( base::mean(X1) - base::mean(X0) )
		xdev = 2
		dev_ok = FALSE
		
		while( !dev_ok )
		{
			dev_ok = TRUE
			xmin = base::min(Y0,X0,X1) - xdev * xdiff
			xmax = base::max(Y0,X0,X1) + xdev * xdiff
			
			x = base::seq( xmin , xmax , length = 200 )
			cdfY0 = rvY0$cdf(x)
			cdfX0 = rvX0$cdf(x)
			cdfX1 = rvX1$cdf(x)
			cdfY1 = rvY0$cdf(rvX0$icdf(cdfX1))
			
			if( base::min(Y0) < base::min(X1s) )
			{
				i = base::max( base::which( x < rvY0$icdf(cdfY1[1]) ) )
				j = base::max( base::which( x < base::min(X1s) ) )
				if( i < j )
				{
					cdfY1[(j-i+1):j] = cdfY0[1:i]
					cdfY1[1:(j-i+1)] = 0.
				}
				else
				{
					cdfY1[1:j] = cdfY0[(i-j+1):i]
				}
			}
			
			if( cdfY1[base::length(cdfY1)] < 1 )
			{
				i = base::max( base::which( x < rvY0$icdf(cdfY1[length(cdfY1)]) ) )
				j = base::min( base::which( cdfY1[1:(length(cdfY1)-1)] == cdfY1[length(cdfY1)] ) )
				if( j > 0 )
				{
					dif = min( length(x) - j , length(x) - i )
					cdfY1[j:(j+dif)] = cdfY0[i:(i+dif)]
					if( j + dif < length(x) ) cdfY1[(j+dif):length(cdfY1)] = 1.
				}
				else
				{
					dev_ok = FALSE
					xdev = 2 * xdev
				}
			}
		}
		
		cdfY1_fct = approxfun( cdfY1 , x )
		Y1 = cdfY1_fct( stats::runif( n = 10000 , min = base::min(cdfY1) , max = base::max(cdfY1) ) )
		
		return(Y1)
	},
	##}}}
	
	infer_Y1_new = function( Y0 , X0 , X1 , idx )##{{{
	{
		dsupp = 1000
		tol   = 1e-6
		samples_Y1 = 10000
		
		## Normalization
		mY0 = base::mean(Y0)
		mX0 = base::mean(X0)
		mX1 = base::mean(X1)
		sY0 = stats::sd(Y0)
		sX0 = stats::sd(X0)
		
		X0s = (X0 - mX0) * sY0 / sX0 + mY0
		X1s = (X1 - mX1) * sY0 / sX0 + mX1 + mY0 - mX0
		
		
		## CDF
		rvY0  = self$distY0$law[[idx]]
		rvX0s = base::do.call( self$distX0$dist[[idx]]$new , self$distX0$kwargs )
		rvX0s$fit(X0s)
		rvX1s = base::do.call( self$distX1$dist[[idx]]$new , self$distX1$kwargs )
		rvX1s$fit(X1s)
		
		## Support
		## Here the support is such that the CDF of Y0, X0s and X1s start from 0
		## and go to 1
		x_min = base::min(Y0,X0s,X1s,X0,X1)
		x_max = base::max(Y0,X0s,X1s,X0,X1)
		x_eps = 0.05 * (x_max - x_min)
		x_fac = 1
		x = base::seq( x_min - x_fac * x_eps , x_max + x_fac * x_eps , length = dsupp )
		
		support_test = function( rv , x )
		{
			if( !base::abs(rv$cdf(x[1])) < tol )
				return(FALSE)
			if( !base::abs(rv$cdf(x[length(x)])-1) < tol )
				return(FALSE)
			return(TRUE)
		}
		
		while( !support_test(rvY0,x) || !support_test(rvX0s,x) || !support_test(rvX1s,x) )
		{
			x_fac = 2 * x_fac
			x = base::seq( x_min - x_fac * x_eps , x_max + x_fac * x_eps , length = dsupp )
		}
		x_fac = x_fac / 2
		dsupp = as.integer(dsupp/1.2)
		
		## Loop to check the support
		extend_support = TRUE
		while(extend_support)
		{
			
			## Support
			extend_support = FALSE
			dsupp = as.integer(dsupp*1.2)
			x_fac = 2*x_fac
			x = base::seq( x_min - x_fac * x_eps , x_max + x_fac * x_eps , length = dsupp )
			
			## Inference of the CDF of Y1
			cdfY1 = rvY0$cdf(rvX0s$icdf(rvX1s$cdf(x)))
			
			## Correction of the CDF, we want that the CDF of Y1 start from 0 and goto 1
			if( !(base::abs(cdfY1[1]) < tol) )
			{
				## CDF not start at 0
				idx  = base::max(which(base::abs(cdfY1[1] - cdfY1) < tol))
				if( idx == 1 )
				{
					extend_support = TRUE
				}
				else
				{
					supp_l_X0s = rvX0s$icdf(cdfY1[1]) - rvX0s$icdf(0)
					supp_l_X1s = rvX1s$icdf(cdfY1[1]) - rvX1s$icdf(0)
					supp_l_Y0  = rvY0$icdf(cdfY1[1])  - rvY0$icdf(0)
					supp_l_Y1  = supp_l_Y0 * supp_l_X1s / supp_l_X0s
					if( x[idx] - supp_l_Y1 < x[1] )
					{
						extend_support = TRUE
					}
					else
					{
						idxl = base::which.min(base::abs(x - (x[idx] - supp_l_Y1)))
						cdfY1[1:(idxl-1)] = 0
						cdfY1[idxl:idx] = rvY0$cdf( base::seq( rvY0$icdf(0) , rvY0$icdf(cdfY1[idx]) , length = idx - idxl + 1 ) )
					}
				}
			}
			
			size = length(cdfY1)
			if( !(abs(cdfY1[size] - 1) < tol) )
			{
				## CDF not finished at 1
				idx = base::min(which(base::abs(cdfY1[size] - cdfY1) < tol))
				if( idx == dsupp )
				{
					extend_support = TRUE
				}
				else
				{
					supp_r_X0s = rvX0s$icdf(1) - rvX0s$icdf(cdfY1[size]) 
					supp_r_X1s = rvX1s$icdf(1) - rvX1s$icdf(cdfY1[size]) 
					supp_r_Y0  = rvY0$icdf(1)  - rvY0$icdf(cdfY1[size])  
					supp_r_Y1  = supp_r_Y0 * supp_r_X1s / supp_r_X0s
					if( x[idx] + supp_r_Y1 > x[size] )
					{
						extend_support = TRUE
					}
					else
					{
						idxr = which.min(base::abs(x - (x[idx] + supp_r_Y1)))
						cdfY1[(idxr+1):size] = 1
						cdfY1[idx:idxr] = rvY0$cdf( base::seq( rvY0$icdf(cdfY1[idx]) , rvY0$icdf(1) , length = idxr - idx + 1 ) )
					}
				}
			}
		}
		
		
		## Draw Y1
		icdfY1 = stats::approxfun( cdfY1 , x , ties = "ordered" )
		hY1    = icdfY1( runif( n = samples_Y1 , min = 0 , max = 1 ) )
		
		return(hY1)
	}
	##}}}
	
	)
	##}}}
)

