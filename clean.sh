#!/bin/sh

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

## Delete python temporary files
rm -rf python/SBCK.egg*
rm -rf python/build
rm -rf python/dist
rm -rf python/tmp
rm -rf python/var


## Delete R temporary files
rm -f R/SBCK/R/RcppExports.R
rm -f R/SBCK/src/*.so
rm -f R/SBCK/src/*.o
rm -f R/SBCK/src/RcppExports.cpp
rm -f R/*.tar.gz
rm -rf R/SBCK.Rcheck
rm -f R/*.pdf
