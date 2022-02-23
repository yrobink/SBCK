
## Copyright(c) 2022 Yoann Robin
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

## PPPFunctionLink ##{{{

#' PPPFunctionLink
#'
#' @description
#' Base class to build link function pre-post processing class
#'
#' @details
#' Base class to build link function pre-post processing class
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#'
#' ## Bias correction
#' ## Step 1 : construction of the class RBC
#' ppp = PPPSSR$new( bc_method = SBCK::CDFt , cols = 2 ) 
#' ## Step 2 : Fit the bias correction model
#' ppp$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = ppp$predict(X1,X0) 
#' ## Z$Z0 # BC of X0
#' ## Z$Z1 # BC of X1
#' @export
PPPFunctionLink = R6::R6Class( "PPPFunctionLink" ,
	
	inherit = PrePostProcessing,
	
	## Public list {{{
	
	public = list(
	
	## initialize ##{{{
	#' @description
    #' Create a new PPPFunctionLink object.
    #' @param transform_ The transform function
    #' @param itransform_ The inverse transform function
    #' @param cols Columns to apply the link function
	#' @param ... Others arguments are passed to PrePostProcessing
    #' @return A new `PPPFunctionLink` object.
	initialize = function( transform_ , itransform_ , cols = NULL , ... )
	{
		base::do.call( super$initialize , list(...) )
		private$.transform  = transform_
		private$.itransform = itransform_
		private$.cols       = cols
	},
	##}}}
	
	## transform ##{{{
	#' @description
    #' Apply the transform.
    #' @param X Data to transform
    #' @return Xt a transformed matrix
	transform = function( X )
	{
		if( is.null(private$.cols) )
			return(private$.transform(X))
		
		Xt = X
		Xt[,private$.cols] = private$.transform(Xt[,private$.cols])
		return(Xt)
	},
	##}}}
	
	## itransform ##{{{
	#' @description
    #' Apply the inverse transform.
    #' @param Xt Data to transform
    #' @return X a transformed matrix
	itransform = function( Xt )
	{
		if( is.null(private$.cols) )
			return(private$.itransform(Xt))
		
		X = Xt
		X[,private$.cols] = private$.itransform(X[,private$.cols])
		return(X)
	}
	##}}}
	
	),
	##}}}
	
	## Private list ##{{{
	
	private = list(
	
	.transform  = NULL,
	.itransform = NULL,
	.cols       = NULL
	)
	##}}}
	
)
##}}}

## PPPSquareLink ##{{{

#' PPPSquareLink
#'
#' @description
#' Square link function
#'
#' @details
#' Square link function
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#'
#' ## Bias correction
#' ## Step 1 : construction of the class RBC
#' ppp = PPPSSR$new( bc_method = SBCK::CDFt , cols = 2 ) 
#' ## Step 2 : Fit the bias correction model
#' ppp$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = ppp$predict(X1,X0) 
#' ## Z$Z0 # BC of X0
#' ## Z$Z1 # BC of X1
#' @export
PPPSquareLink = R6::R6Class( "PPPSquareLink" ,
	
	inherit = PPPFunctionLink,
	
	## Public list {{{
	
	public = list(
	
	## initialize ##{{{
	#' @description
    #' Create a new PPPSquareLink object.
    #' @param cols Columns to apply the link function
	#' @param ... Others arguments are passed to PrePostProcessing
    #' @return A new `PPPSquareLink` object.
	initialize = function( cols = NULL , ... )
	{
		kwargs = list(...)
		kwargs[["transform_"]]  = function(x) { return(x^2) }
		kwargs[["itransform_"]] = function(x) { return(where( x > 0 , base::sqrt(base::abs(x)) , -base::sqrt(base::abs(x)) ) ) }
		kwargs[["cols"]]        = cols
		base::do.call( super$initialize , kwargs )
	}
	##}}}
	
	),
	##}}}
	
	## Private list ##{{{
	
	private = list(
	
	)
	##}}}
	
)
##}}}

## PPPLogLinLink  ##{{{

#' PPPLogLinLink
#'
#' @description
#' Log linear link function
#'
#' @details
#' Log linear link function
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#'
#' ## Bias correction
#' ## Step 1 : construction of the class RBC
#' ppp = PPPSSR$new( bc_method = SBCK::CDFt , cols = 2 ) 
#' ## Step 2 : Fit the bias correction model
#' ppp$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = ppp$predict(X1,X0) 
#' ## Z$Z0 # BC of X0
#' ## Z$Z1 # BC of X1
#' @export
PPPLogLinLink = R6::R6Class( "PPPLogLinLink" ,
	
	inherit = PPPFunctionLink,
	
	## Public list {{{
	
	public = list(
	
	## initialize ##{{{
	#' @description
    #' Create a new PPPLogLinLink object.
    #' @param cols Columns to apply the link function
	#' @param ... Others arguments are passed to PrePostProcessing
    #' @return A new `PPPLogLinLink` object.
	initialize = function( cols = NULL , ... )
	{
		kwargs = list(...)
		kwargs[["transform_"]]  = function(x) { return( where( (0 < x) & (x < 1) , base::log(x) , x - 1 ) ) }
		kwargs[["itransform_"]] = function(x) { return( where( x < 0 , base::exp(x) , x + 1 ) ) }
		kwargs[["cols"]]        = cols
		base::do.call( super$initialize , kwargs )
	}
	##}}}
	
	),
	##}}}
	
)
##}}}

