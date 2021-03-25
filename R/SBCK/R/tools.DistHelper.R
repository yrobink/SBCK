
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


###############
## Libraries ##
###############

###############
## Functions ##
###############


#' Dist Helper
#'
#' @description
#' Class used by CDFt and QM to facilitate fit, do not use
#'
#' @details
#' Used to parallel work for margins
#'
#' @examples
#' ## 
#' @export
DistHelper = R6::R6Class( "DistHelper" ,
	
	public = list(
	
	#' @field dist [ROOPSD distribution] name of class
	dist   = NULL,
	#' @field law [ROOPSD distribution] class set
	law    = NULL,
	#' @field kwargs [list] arguments of dist
	kwargs = NULL,
	
	## initialize ##{{{
	#' @description
    #' Create a new DistHelper object.
	#' @param dist [ROOPSD distribution or list] statistical law
	#' @param kwargs [list] arguments passed to dist
	#'
	#' @return A new `DistHelper` object.
	initialize = function( dist , kwargs )
	{
		self$dist   = if( is.null(dist) ) ROOPSD::rv_histogram else dist
		self$kwargs = if( is.null(kwargs) ) list() else kwargs
		self$law    = list()
	},
	##}}}
	
	## set_features ##{{{
	#' @description
    #' set the number of features
    #' @param n_features [integer] numbers of features
	#'
    #' @return NULL
	set_features = function( n_features )
	{
		if( !is.list(self$dist) )
		{
			dist = list()
			for( i in 1:n_features )
				dist[[i]] = self$dist
			self$dist = dist
		}
	},
	##}}}
	
	## fit ##{{{
	#' @description
	#' fit the laws
	#' @param X [matrix] dataset to fit
	#' @param i [integer] margins to fit
	#'
    #' @return NULL
	fit = function( X , i )
	{
		if( "R6" %in% class(self$dist[[i]]) )
		{
			self$law[[i]] = self$dist[[i]]
		}
		else
		{
			self$law[[i]] = base::do.call( self$dist[[i]]$new , self$kwargs )
			self$law[[i]]$fit(X)
		}
	},
	##}}}
	
	## is_frozen ##{{{
	#' @description
	#' Test if margins i is frozen
	#' @param i [integer] margins to fit
	#'
	#' @return [bool]
	is_frozen = function(i)
	{
		return( "R6" %in% class(self$dist[[i]]) )
	},
	##}}}
	
	## is_parametric ##{{{
	#' @description
	#' Test if margins i is parametric
	#' @param i [integer] margins to fit
	#'
	#' @return [bool]
	is_parametric = function(i)
	{
		return( "AbstractDist" %in% class(self$law[[i]]) )
	}
	## }}}
	
	)
)




