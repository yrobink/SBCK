
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

#' where function
#'
#' This function return a vector / matrix / array of the same shape than
#' cond / x / y such that if(cond) values are x, and else y.
#' 
#' @usage where(cond,x,y)
#' @param cond [vector/matrix/array] Boolean values
#' @param x [vector/matrix/array] Values if cond is TRUE
#' @param y [vector/matrix/array] Values if cond is FALSE
#'
#' @return z [vector/matrix/array].
#'
#' @examples
#' x = base::seq( -2 , 2 , length = 100 )
#' y = where( x < 1 , x , exp(x) ) ## y = x if x < 1, else exp(x)
#'
#' @export
where = function( cond , x , y )
{
	if( length(y) == 1 )
	{
		y = x * 0 + y
	}
	if( length(x) == 1 )
	{
		x = y * 0 + x
	}
	
	z = x
	z[!cond] = y[!cond]
	return(z)
}

