
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


//=========//
// Include //
//=========//

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "SparseHist.hpp"
#include "NetworkSimplex.hpp"

//============//
// namespaces //
//============//

namespace py = pybind11 ;


//========//
// Module //
//========//

PYBIND11_MODULE( __tools_cpp , m )
{
	//===========//
	// Functions //
	//===========//
	
	m.def( "network_simplex" , &network_simplex ) ;
	
	//=======//
	// Class //
	//=======//
	
	py::class_<SparseHist>( m , "SparseHist" )
	.def( py::init<Eigen::Ref<const Eigen::MatrixXd>, Eigen::Ref<const Eigen::VectorXd>>()                                   , py::arg("X") , py::arg("bin_width") )
	.def( py::init<Eigen::Ref<const Eigen::MatrixXd>, py::list>()                                                            , py::arg("X") , py::arg("bin_width") )
	.def( py::init<Eigen::Ref<const Eigen::MatrixXd>, Eigen::Ref<const Eigen::VectorXd>,Eigen::Ref<const Eigen::VectorXd>>() , py::arg("X") , py::arg("bin_width") , py::arg("bin_origin") )
	.def( py::init<Eigen::Ref<const Eigen::MatrixXd>, py::list,Eigen::Ref<const Eigen::VectorXd>>()                          , py::arg("X") , py::arg("bin_width") , py::arg("bin_origin") )
	.def( py::init<Eigen::Ref<const Eigen::MatrixXd>, Eigen::Ref<const Eigen::VectorXd>,py::list>()                          , py::arg("X") , py::arg("bin_width") , py::arg("bin_origin") )
	.def( py::init<Eigen::Ref<const Eigen::MatrixXd>, py::list,py::list>()                                                   , py::arg("X") , py::arg("bin_width") , py::arg("bin_origin") )
	.def( "__repr__" , &SparseHist::repr )
	.def( "argwhere" , &SparseHist::argwhere , py::arg("X") )
	.def_readwrite( "bin_width"  , &SparseHist::m_bin_width   )
	.def_readwrite( "bin_origin" , &SparseHist::m_bin_origin  )
	.def_readwrite( "p"          , &SparseHist::m_p           )
	.def_readwrite( "c"          , &SparseHist::m_c           )
	.def_readwrite( "dim"        , &SparseHist::m_dim         )
	.def_readwrite( "size"       , &SparseHist::m_size        )
	;
//	.def_property( "name" , &Pet::getName , &Pet::setName ) ;
	
	
	//============//
	// Attributes //
	//============//
	
//	m.doc() = "pybind11 example plugin" ; // optional module docstring
	m.attr("__name__") = "SBCK.tools.__tools_cpp";
//	m.attr("the_answer") = 42;
}
