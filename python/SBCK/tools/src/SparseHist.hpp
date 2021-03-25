
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

#ifndef SBCK_SPARSEHIST
#define SBCK_SPARSEHIST

//-----------//
// Libraries //
//-----------//

#include <iostream>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <cmath>
#include <pybind11/pybind11.h>


//============//
// namespaces //
//============//

namespace py = pybind11 ;


//=======//
// Class //
//=======//

struct SparseHist
{
	// Some typedef {{{
	typedef unsigned int size_type ;
	typedef Eigen::VectorXi VectIndex ;
	typedef Eigen::VectorXd VectValue ;
	typedef Eigen::MatrixXd DataType ;
	typedef std::function<bool(const VectIndex&,const VectIndex&)> OrderType ;
	typedef std::map<VectIndex,int,OrderType> HashTable ;
	
	//}}}
	
	// Constructor / Destructor {{{
	
	SparseHist( Eigen::Ref<const DataType> X , Eigen::Ref<const VectValue> bin_width ): //{{{
		m_dim(bin_width.size()) ,
		m_size(0) ,
		m_bin_width(bin_width) ,
		m_bin_origin( Eigen::ArrayXd::Zero(m_dim) ) ,
		m_alpha() ,
		m_beta() ,
		m_map( []( const VectIndex& x , const VectIndex& y ) { return std::lexicographical_compare( &x(0) , &x(0) + x.size() , &y(0) , &y(0) + y.size() ) ; } ) ,
//		m_map( []( const VectIndex& x , const VectIndex& y ) { return std::lexicographical_compare( x.begin() , x.end() , y.begin() , y.end() ) ; } ) ,
		m_c() ,
		m_p()
	{
		initialize( X ) ;
	}
	//}}}
	
	SparseHist( Eigen::Ref<const DataType> X , py::list bin_width ): //{{{
		m_dim(py::len(bin_width)) ,
		m_size(0) ,
		m_bin_width(m_dim) ,
		m_bin_origin( Eigen::ArrayXd::Zero(m_dim) ) ,
		m_alpha() ,
		m_beta() ,
		m_map( []( const VectIndex& x , const VectIndex& y ) { return std::lexicographical_compare( &x(0) , &x(0) + x.size() , &y(0) , &y(0) + y.size() ) ; } ) ,
		m_c() ,
		m_p()
	{
		int s(0) ;
		for( auto item : bin_width )
			m_bin_width[s++] = item.cast<double>() ;
		
		initialize( X ) ;
	}
	//}}}
	
	SparseHist( Eigen::Ref<const DataType> X , Eigen::Ref<const VectValue> bin_width , Eigen::Ref<const VectValue> bin_origin )://{{{
		m_dim(bin_width.size()) ,
		m_size(0) ,
		m_bin_width(bin_width) ,
		m_bin_origin(bin_origin) ,
		m_alpha() ,
		m_beta() ,
		m_map( []( const VectIndex& x , const VectIndex& y ) { return std::lexicographical_compare( &x(0) , &x(0) + x.size() , &y(0) , &y(0) + y.size() ) ; } ) ,
		m_c() ,
		m_p()
	{
		initialize( X ) ;
	}
	//}}}
	
	SparseHist( Eigen::Ref<const DataType> X , py::list bin_width , py::list bin_origin )://{{{
		m_dim(py::len(bin_width)) ,
		m_size(0) ,
		m_bin_width(m_dim) ,
		m_bin_origin(m_dim) ,
		m_alpha() ,
		m_beta() ,
		m_map( []( const VectIndex& x , const VectIndex& y ) { return std::lexicographical_compare( &x(0) , &x(0) + x.size() , &y(0) , &y(0) + y.size() ) ; } ) ,
		m_c() ,
		m_p()
	{
		int s(0) ;
		for( auto item : bin_width )
			m_bin_width[s++] = item.cast<double>() ;
		
		s = 0 ;
		for( auto item : bin_origin )
			m_bin_origin[s++] = item.cast<double>() ;
		
		initialize( X ) ;
	}
	//}}}
	
	SparseHist( Eigen::Ref<const DataType> X , py::list bin_width , Eigen::Ref<const VectValue> bin_origin )://{{{
		m_dim(bin_origin.size()) ,
		m_size(0) ,
		m_bin_width(m_dim) ,
		m_bin_origin(bin_origin) ,
		m_alpha() ,
		m_beta() ,
		m_map( []( const VectIndex& x , const VectIndex& y ) { return std::lexicographical_compare( &x(0) , &x(0) + x.size() , &y(0) , &y(0) + y.size() ) ; } ) ,
		m_c() ,
		m_p()
	{
		int s(0) ;
		for( auto item : bin_width )
			m_bin_width[s++] = item.cast<double>() ;
		
		initialize( X ) ;
	}
	//}}}
	
	SparseHist( Eigen::Ref<const DataType> X , Eigen::Ref<const VectValue> bin_width , py::list bin_origin )://{{{
		m_dim(bin_width.size()) ,
		m_size(0) ,
		m_bin_width(bin_width) ,
		m_bin_origin(m_dim) ,
		m_alpha() ,
		m_beta() ,
		m_map( []( const VectIndex& x , const VectIndex& y ) { return std::lexicographical_compare( &x(0) , &x(0) + x.size() , &y(0) , &y(0) + y.size() ) ; } ) ,
		m_c() ,
		m_p()
	{
		int s(0) ;
		for( auto item : bin_origin )
			m_bin_origin[s++] = item.cast<double>() ;
		
		initialize( X ) ;
	}
	//}}}
	
	void initialize( const DataType& X ) //{{{
	{
		// Linear mapping
		m_alpha      = 1. / m_bin_width.array() ;
		m_beta       = - m_bin_origin.array() * m_alpha.array() ;
		
		// Bins estimation
		for( int s = 0 ; s < X.rows() ; ++s )
		{
			m_map[bin_index(X.row(s))]++ ;
		}
		
		
		// Final construction
		m_size = m_map.size() ;
		double dsize = static_cast<double>(X.rows()) ;
		m_p.resize( m_size ) ;
		m_c.resize( m_size , m_dim ) ;
		size_type s = 0 ;
		for( auto& keyval : m_map )
		{
			m_p[s] = keyval.second / dsize ;
			m_c.row(s++) = bin_center(keyval.first) ;
		}
	}
	//}}}
	
	~SparseHist()//{{{
	{}
	//}}}
	
	// }}}
	
	std::string repr()//{{{
	{
		std::string _repr("") ;
		_repr += "SBCK.tools.SparseHist\n" ;
		_repr += "=====================\n" ;
		_repr += "* size:" + std::to_string(m_size) + "\n" ;
		_repr += "* dim :" + std::to_string(m_dim)  + "\n" ;
		return _repr ;
	}
	//}}}
	
	// Methods {{{
	
	VectIndex bin_index( const VectValue& x )
	{
		VectIndex index = ( m_alpha.array() * x.array() + m_beta.array() ).floor().cast<int>() ;
		return index ;
	}
	
	VectValue bin_center( const VectIndex& index )
	{
		VectValue x(m_dim) ;
		for( size_type s = 0 ; s < m_dim ; ++s )
		{
			x[s] = m_bin_origin[s] + m_bin_width[s] * static_cast<double>(index[s]) + m_bin_width[s] / 2. ;
		}
		return x ;
	}
	
	VectIndex argwhere( Eigen::Ref<const DataType> X )
	{
		VectIndex index ;
		VectIndex lIndex(Eigen::VectorXi::Zero(X.rows())) ;
		typename HashTable::iterator it ;
		for( int s = 0 ; s < X.rows() ; ++s )
		{
			index = bin_index(X.row(s)) ;
			it = m_map.find( index ) ;
			lIndex[s] = static_cast<int>( ( it == m_map.end() ) ? -1 : std::distance( m_map.begin() , it ) ) ;
		}
		return lIndex ;
	}
	//}}}
	
	// Arguments {{{
	size_type	m_dim ;
	size_type	m_size ;
	VectValue	m_bin_width ;
	VectValue	m_bin_origin ;
	VectValue	m_alpha ;
	VectValue	m_beta ;
	HashTable	m_map ;
	DataType	m_c ;
	VectValue	m_p ;
	//}}}
};


#endif

