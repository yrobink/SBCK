
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

## PPPSSR ##{{{

#' PPPSSR
#'
#' @description
#' Apply the SSR transformation.
#'
#' @details
#' Apply the SSR transformation. The SSR transformation replace the 0 by a
#' random values between 0 and the minimal non zero value (the threshold). The
#' inverse transform replace all values lower than the threshold by 0. The
#' threshold used for inverse transform is given by the keyword `isaved`, which
#' takes the value `Y0` (reference in calibration period), or `X0` (biased in
#' calibration period), or `X1` (biased in projection period)
#'
#' @examples
#' ## Start with data
#' XY = SBCK::dataset_like_tas_pr(2000)
#' X0 = XY$X0
#' X1 = XY$X1
#' Y0 = XY$Y0
#' 
#' ## Define the PPP method
#' ppp = PPPSSR$new( bc_method = CDFt , cols = 2 )
#' 
#' ## And now the correction
#' ## Bias correction
#' ppp$fit(Y0,X0,X1)
#' Z = ppp$predict(X1,X0)
#' 
#' @export
PPPSSR = R6::R6Class( "PPPSSR" ,
	
	inherit = PrePostProcessing,
	
	## Public list {{{
	
	public = list(
	
	## Public arguments
	#' @field Xn [vector] Threshold
	Xn = NULL,
	
	## initialize ##{{{
	#' @description
    #' Create a new PPPSSR object.
    #' @param cols Columns to apply the SSR
	#' @param isaved Choose the threshold used for the inverse transform. Can be
	#'        "Y0", "X0" and "X1".
	#' @param ... Others arguments are passed to PrePostProcessing
    #' @return A new `PPPSSR` object.
	initialize = function( cols = NULL , isaved = "Y0" , ...)
	{
		kwargs = list(...)
		base::do.call( super$initialize , kwargs )
		private$.cols     = cols
		private$.isaved   = isaved
		private$.icurrent = 0
		
		if( isaved == "Y0" )
			private$.isaved = 1
		if( isaved == "X0" )
			private$.isaved = 2
		if( isaved == "X1" )
			private$.isaved = 3
	},
	##}}}
	
	## transform  ##{{{
	#' @description
    #' Apply the SSR transform, i.e. all 0 are replaced by random values between
    #' 0 (excluded) and the minimal non zero value.
    #' @param X Data to transform
    #' @return Xt a transformed matrix
	transform = function( X )
	{
		private$.icurrent = private$.icurrent + 1
		Xt = X
		if( is.null(private$.cols) )
		{
			private$.cols = 1:ncol(X)
		}
		cols = private$.cols
		
		Xn = base::apply( where( Xt[,cols,drop=FALSE] > 0 , Xt[,cols,drop=FALSE] , NaN ) , 2 , base::min , na.rm = TRUE )
		if( private$.icurrent == private$.isaved )
			self$Xn = Xn
		
		ncols = length(cols)
		A = where( Xt[,cols,drop=FALSE] > 0 , Xt[,cols,drop=FALSE] , base::t( matrix( stats::runif( nrow(Xt) * ncols , min = Xn / 1000 , max = Xn ) , nrow = ncols ) ) )
		Xt[,cols] = where( Xt[,cols,drop=FALSE] > 0 , Xt[,cols,drop=FALSE] , base::t( matrix( stats::runif( nrow(Xt) * ncols , min = Xn / 1000 , max = Xn ) , nrow = ncols ) ) )
		return(Xt)
	},
	##}}}
	
	## itransform  ##{{{
	#' @description
    #' Apply the inverse SSR transform, i.e. all values lower than the threshold
    #' found in the transform function are replaced by 0.
    #' @param Xt Data to transform
    #' @return X a transformed matrix
	itransform = function( Xt )
	{
		cols = private$.cols
		X = Xt
		X[,cols] = where( base::t( base::t(X[,cols,drop=FALSE]) < self$Xn ) , 0 , X[,cols,drop=FALSE] )
		return(X)
	}
	##}}}
	
	),
	##}}}
	
	## Private list ##{{{
	
	private = list(
	
	## Private arguments ##{{{
	.cols     = NULL,
	.isaved   = NULL,
	.icurrent = NULL
	
	##}}}
	
	## Private methods ##{{{
	
	
	##}}}
	
	)
	##}}}
	
)
##}}}

