 
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


#ifndef SBCK_NETWORKSIMPLEX
#define SBCK_NETWORKSIMPLEX



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

#include "NetworkSimplexLemon.hpp"


//============//
// namespaces //
//============//

namespace py = pybind11 ;


//===========//
// Functions //
//===========//

using namespace lemon ;
typedef unsigned int node_id_type ;
typedef Eigen::VectorXd VectValue ;
typedef Eigen::MatrixXd DataType ;


//void network_simplex( int size1 , int size2 , double* X , double* Y , double* D , double* G , double* cost )
py::tuple network_simplex( Eigen::Ref<const VectValue> X , Eigen::Ref<const VectValue> Y , Eigen::Ref<const DataType> D )
{
	// beware M and C anre strored in row major C style!!!
	int cur ;
	double max , max_iter(-1.) ;
	node_id_type n , m , n1(X.size()) , n2(Y.size()) ;
	
	typedef FullBipartiteDigraph Digraph ;
	DIGRAPH_TYPEDEFS(FullBipartiteDigraph) ;
	
	// Get the number of non zero coordinates for r and c
	n = 0 ;
	for( node_id_type i = 0 ; i < n1 ; i++ )
	{
		if( X(i) > 0 )
			n++ ;
	}
	m = 0 ;
	for( node_id_type i = 0 ; i < n2 ; i++ )
	{
		if( Y(i) > 0 )
			m++ ;
	}
	
	Eigen::MatrixXd G(n1,n2) ;
	
	// Define the graph
	std::vector<int> indI(n) , indJ(m) ;
	std::vector<double> weights1(n) , weights2(m) ;
	Digraph di(n,m) ;
	NetworkSimplexSimple<Digraph,double,double,node_id_type> net( di , true , n + m , n * m , max_iter ) ;
	

	// Set supply and demand, don't account for 0 values (faster)
	max = 0 ;
	cur = 0 ;
	for( node_id_type i = 0 ; i < n1 ; i++ )
	{
		double val = X(i) ;
		if( val > 0 )
		{
			weights1[ di.nodeFromId(cur) ] = val ;
			max += val ;
			indI[cur++] = i ;
		}
	}
	

	// Demand is actually negative supply...
	max = 0 ;
	cur = 0 ;
	for( node_id_type i = 0 ; i < n2 ; i++ )
	{
		double val = Y(i) ;
		if( val > 0 )
		{
			weights2[ di.nodeFromId(cur) ] = -val ;
			indJ[cur++] = i ;
			max -= val ;
		}
	}
	
	
	net.supplyMap( &weights1[0] , n , &weights2[0] , m ) ;
	

	// Set the cost of each edge
	max = 0 ;
	for( node_id_type i = 0 ; i < n ; i++ )
	{
		for( node_id_type j = 0 ; j < m ; j++ )
		{
			double val = D(indI[i],indJ[j]) ;//*( D + indI[i] * n2 + indJ[j] ) ;
			net.setCost( di.arcFromId( i * m + j ) , val ) ;
			if( val > max )
			{
				max = val ;
			}
		}
	}
	
	
	// Solve the problem with the network simplex algorithm
	int ret = net.run() ;
	if( ret == static_cast<int>(net.OPTIMAL) )
	{
		for( node_id_type i = 0 ; i < n ; i++ )
		{
			for( node_id_type j = 0 ; j < m ; j++ )
			{
				G(indI[i],indJ[j]) = net.flow( di.arcFromId( i * m + j ) ) ;
			}
		}
	}
	
	py::tuple output = py::make_tuple( G , ret ) ;
	
	return output ;
}


#endif
