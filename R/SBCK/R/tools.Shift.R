
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

#' Shift
#'
#' @description
#' Class to shift a dataset.
#'
#' @details
#' Transform autocorrelations to intervariables correlations
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(lag,method,ref,)}}{This method is used to create object of this class with \code{Shift}}
#'   \item{\code{transform(X)}}{Method to shift a dataset}
#'   \item{\code{inverse(Xs)}}{Method to inverse the shift of a dataset}
#' }
#' @examples
#' X = base::t(matrix( 1:20 , nrow = 2 , ncol = 10 ))
#' 
#' sh = Shift$new(1)
#' Xs = sh$transform(X)
#' Xi = sh$inverse(Xs)
#'
#' @export
Shift = R6::R6Class( "Shift" ,
	
	private = list( ##{{{
	
	.method = NULL,
	.ref    = NULL,
	
	inverse_by_row = function(Xs)##{{{
	{
		n_features = as.integer( dim(Xs)[2] / ( self$lag + 1 ))
		n_samples  = dim(Xs)[1] + self$lag
		Xi = matrix( NA , ncol = n_features , nrow = n_samples )
		for( r in base::c(1:self$lag,self$ref) )
		{
			idx = base::seq( r , n_samples - self$lag , self$lag )
			Xs0 = as.vector(base::t(Xs[idx[-length(idx)],1:(base::ncol(Xs) - n_features)]))
			Xs1 = as.vector(base::t(Xs[idx[length(idx)],]))
			Xs01 = base::t(matrix( base::c(Xs0,Xs1) , nrow = n_features ))
			Xi[r:(r+base::nrow(Xs01)-1),] = Xs01
		}
		
		return(Xi)
	},
	##}}}
	
	inverse_by_col = function(Xs)##{{{
	{
		n_features = as.integer( dim(Xs)[2] / (self$lag + 1) )
		n_samples  = dim(Xs)[1] + self$lag
		Xi = matrix( NA , nrow = n_samples , ncol = n_features )
		
		ref = self$ref
		if( ref < 1 )
			ref = 1
		if( ref > n_features )
			ref = n_features
		
		for( i in base::c( 0:self$lag , ref ) )
		{
			db = i * n_features + 1
			de = i * n_features + n_features
			tb = i + 1
			te = n_samples - ( self$lag + 1 ) + i + 1
			Xi[tb:te,] = Xs[,db:de]
		}
		
		return(Xi)
	}
	##}}}
	
	),
	##}}}
	
	active = list( ##{{{
	
	method = function(value)##{{{
	{
		if( missing(value) )
		{
			return(private$.method)
		}
		else
		{
			if( value %in% base::c("col","row") )
				private$.method = value
			else
				private$.method = "row"
		}
	},
	##}}}
	
	ref = function(value)##{{{
	{
		if( missing(value) )
		{
			return(private$.ref)
		}
		else
		{
			private$.ref = ( (value - 1) %% ( self$lag + 1 ) ) + 1
		}
	}
	
	),
	##}}}
	
	##}}}
	
	public = list( ##{{{
	
	###############
	## Arguments ##
	###############
	
	#' @field lag [integer] max lag for autocorrelations
	lag    = NULL,
	#' @field method [character] If inverse is by row or column.
	#' @field ref [integer] reference column/row to inverse shift.
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new Shift object.
	#'
	#' @param lag [integer] max lag for autocorrelations
	#' @param method [character] If "row" inverse by row, else by column
	#' @param ref [integer] starting point for inverse transform
	#'
	#' @return A new `Shift` object.
	initialize = function( lag , method = "row" , ref = 1 )
	{
		self$lag    = lag
		self$ref    = ref
		self$method = method
	},
	##}}}
	
	
	#############
	## Methods ##
	#############
	
	## transform ##{{{
	#' @description
    #' Shift the data
    #' @param X [matrix: n_samples * n_features] Data to shift
    #'
    #' @return [matrix] Matrix shifted
	transform = function( X )
	{
		if( is.vector(X) )
			X = matrix( X , nrow = length(X) , ncol = 1 )
		
		n_samples  = base::dim(X)[1]
		n_features = base::dim(X)[2]
		Xs = matrix( NA , nrow = n_samples - self$lag , ncol = ( self$lag + 1 ) * n_features )
		
		for( i in 0:self$lag )
		{
			db = i * n_features + 1
			de = i * n_features + n_features
			tb = i + 1
			te = n_samples - ( self$lag + 1 ) + i + 1
			Xs[,db:de] = X[tb:te,]
		}
		
		return(Xs)
	},
	##}}}
	
	## inverse ##{{{
	#' @description
    #' Inverse the shift of the data
    #' @param Xs [matrix] Data Shifted
    #'
    #' @return [matrix] Matrix un shifted
	inverse = function( Xs )
	{
		if( self$method == "row" )
			return(private$inverse_by_row(Xs))
		else
			return(private$inverse_by_col(Xs))
	}
	##}}}
	
	)
	##}}}
)

