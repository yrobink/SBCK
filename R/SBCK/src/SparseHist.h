
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
#include <Eigen/Dense>


//============//
// namespaces //
//============//



//=======//
// Class //
//=======//

//' @export SparseHistBase
class SparseHistBase
{
	public:
	// Some typedef {{{
	typedef unsigned int size_type ;
	typedef Eigen::VectorXi VectIndex ;
	typedef Eigen::VectorXd VectValue ;
	typedef Eigen::MatrixXd DataType ;
	typedef std::function<bool(VectIndex,VectIndex)> OrderType ;
	typedef std::map<VectIndex,int,OrderType> HashTable ;
	
	//}}}
	
	// Constructor / Destructor {{{
	
	SparseHistBase( DataType X , VectValue bin_width , VectValue bin_origin )://{{{
		m_n_features(bin_width.size()) ,
		m_n_samples(0) ,
		m_bin_width(bin_width) ,
		m_bin_origin(bin_origin) ,
		m_alpha() ,
		m_beta() ,
		m_map( []( VectIndex x , VectIndex y ) { return std::lexicographical_compare( &x(0) , &x(0) + x.size() , &y(0) , &y(0) + y.size() ) ; } ) ,
		m_c() ,
		m_p()
	{
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
		m_n_samples = m_map.size() ;
		double dsize = static_cast<double>(X.rows()) ;
		m_p.resize( m_n_samples ) ;
		m_c.resize( m_n_samples , m_n_features ) ;
		size_type s = 0 ;
		for( auto& keyval : m_map )
		{
			m_p[s] = keyval.second / dsize ;
			m_c.row(s++) = bin_center(keyval.first) ;
		}
	}
	//}}}
	
	~SparseHistBase()//{{{
	{}
	//}}}
	
	// }}}
	
	// Methods {{{
	
	VectIndex bin_index( const VectValue& x )
	{
		VectIndex index = ( m_alpha.array() * x.array() + m_beta.array() ).floor().cast<int>() ;
		return index ;
	}
	
	VectValue bin_center( const VectIndex& index )
	{
		VectValue x(m_n_features) ;
		for( size_type s = 0 ; s < m_n_features ; ++s )
		{
			x[s] = m_bin_origin[s] + m_bin_width[s] * static_cast<double>(index[s]) + m_bin_width[s] / 2. ;
		}
		return x ;
	}
	
	VectIndex argwhere( DataType X )
	{
		VectIndex index ;
		VectIndex lIndex(Eigen::VectorXi::Zero(X.rows())) ;
		typename HashTable::iterator it ;
		for( int s = 0 ; s < X.rows() ; ++s )
		{
			index = bin_index(X.row(s)) ;
			it = m_map.find( index ) ;
			lIndex[s] = static_cast<int>( ( it == m_map.end() ) ? -1 : std::distance( m_map.begin() , it ) ) + 1 ;
		}
		return lIndex ;
	}
	//}}}
	
	
	// Arguments {{{
	size_type	m_n_features ;
	size_type	m_n_samples ;
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
