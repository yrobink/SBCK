// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Copyright(c) 2021 Yoann Robin
// 
// This file is part of SBCK.
// 
// SBCK is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// SBCK is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with SBCK.  If not, see <https://www.gnu.org/licenses/>.


#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
#include "SparseHist.h"
#include "NetworkSimplex.h"
#include "Pairwise_distances.h"


RCPP_MODULE(SBCK_cpp){
	
	// Classes
	Rcpp::class_<SparseHistBase>("SparseHistBase")
	.constructor<Eigen::MatrixXd,Eigen::VectorXd,Eigen::VectorXd>()
	.field( "c"          , &SparseHistBase::m_c          , "" )
	.field( "p"          , &SparseHistBase::m_p          , "" )
	.field( "bin_width"  , &SparseHistBase::m_bin_width  , "" )
	.field( "bin_origin" , &SparseHistBase::m_bin_origin , "" )
	.field( "n_features" , &SparseHistBase::m_n_features , "" )
	.field( "n_samples"  , &SparseHistBase::m_n_samples  , "" )
	.method( "argwhere"  , &SparseHistBase::argwhere     , "" )
	;
}



