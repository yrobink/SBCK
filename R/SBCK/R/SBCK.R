
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

#' SBCK
#'
#' Statistical Bias Correction Kit
#'
#' @docType package
#' @author Yoann Robin Maintainer: Yoann Robin <yoann.robin.k@gmail.com>
#' @import methods
#' @import ROOPSD
#' @import R6
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib SBCK, .registration=TRUE
#' @name SBCK
NULL


loadModule("SBCK_cpp", TRUE)


