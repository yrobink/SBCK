
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
#' Base class to build link function pre-post processing class. See also
#' the PrePostProcessing documentation
#'
#' @details
#' This class is used to define pre/post processing class with a link function
#' and its inverse. See example.
#'
#' @examples
#' ## Start with data
#' XY = SBCK::dataset_like_tas_pr(2000)
#' X0 = XY$X0
#' X1 = XY$X1
#' Y0 = XY$Y0
#' 
#' ## Define the link function
#' transform  = function(x) { return(x^3) }
#' itransform = function(x) { return(x^(1/3)) }
#' 
#' ## And the PPP method
#' ppp = PPPFunctionLink$new( bc_method = CDFt , transform = transform ,
#'                                              itransform = itransform )
#' 
#' ## And now the correction
#' ## Bias correction
#' ppp$fit(Y0,X0,X1)
#' Z = ppp$predict(X1,X0)
#' 
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
#' Square link function. See also the PrePostProcessing documentation.
#'
#' @details
#' Square link function. The transform is x^2, and the sign(x)*sqrt(abs(x)) its
#' inverse.
#'
#' @examples
#' ## Start with data
#' XY = SBCK::dataset_like_tas_pr(2000)
#' X0 = XY$X0
#' X1 = XY$X1
#' Y0 = XY$Y0
#' 
#' ## Define the PPP method
#' ppp = PPPSquareLink$new( bc_method = CDFt , cols = 2 )
#' 
#' ## And now the correction
#' ## Bias correction
#' ppp$fit(Y0,X0,X1)
#' Z = ppp$predict(X1,X0)
#' 
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
#' Log linear link function. See also the PrePostProcessing documentation.
#'
#' @details
#' Log linear link function. The transform is log(x) if 0 < x < 1, else x -1,
#' and the inverse transform exp(x) if x < 0, else x + 1.
#'
#' @examples
#' ## Start with data
#' XY = SBCK::dataset_like_tas_pr(2000)
#' X0 = XY$X0
#' X1 = XY$X1
#' Y0 = XY$Y0
#' 
#' ## Define the PPP method
#' ppp = PPPLogLinLink$new( bc_method = CDFt , cols = 2 ,
#'                          pipe = list(PPPSSR),
#'                          pipe_kwargs = list(list(cols=2)) )
#' 
#' ## And now the correction
#' ## Bias correction
#' ppp$fit(Y0,X0,X1)
#' Z = ppp$predict(X1,X0)
#' 
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

