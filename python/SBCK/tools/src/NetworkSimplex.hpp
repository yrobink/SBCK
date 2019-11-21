 
//==============================================================================//
//==============================================================================//
//                                                                              //
// Copyright Yoann Robin, 2019                                                  //
//                                                                              //
// yoann.robin.k@gmail.com                                                      //
//                                                                              //
// This software is a computer program that is part of the SBCK (Statistical    //
// Bias Correction Kit). This library makes it possible to perform bias         //
// correction with non parametric methods, and give some metrics between Sparse //
// Histogram is high dimensions.                                                //
//                                                                              //
// This software is governed by the CeCILL-C license under French law and       //
// abiding by the rules of distribution of free software.  You can  use,        //
// modify and/ or redistribute the software under the terms of the CeCILL-C     //
// license as circulated by CEA, CNRS and INRIA at the following URL            //
// "http://www.cecill.info".                                                    //
//                                                                              //
// As a counterpart to the access to the source code and  rights to copy,       //
// modify and redistribute granted by the license, users are provided only      //
// with a limited warranty  and the software's author,  the holder of the       //
// economic rights,  and the successive licensors  have only  limited           //
// liability.                                                                   //
//                                                                              //
// In this respect, the user's attention is drawn to the risks associated       //
// with loading,  using,  modifying and/or developing or reproducing the        //
// software by the user in light of its specific status of free software,       //
// that may mean  that it is complicated to manipulate,  and  that  also        //
// therefore means  that it is reserved for developers  and  experienced        //
// professionals having in-depth computer knowledge. Users are therefore        //
// encouraged to load and test the software's suitability as regards their      //
// requirements in conditions enabling the security of their systems and/or     //
// data to be ensured and,  more generally, to use and operate it in the        //
// same conditions as regards security.                                         //
//                                                                              //
// The fact that you are presently reading this means that you have had         //
// knowledge of the CeCILL-C license and that you accept its terms.             //
//                                                                              //
//==============================================================================//
//==============================================================================//

//==============================================================================//
//==============================================================================//
//                                                                              //
// Copyright Yoann Robin, 2019                                                  //
//                                                                              //
// yoann.robin.k@gmail.com                                                      //
//                                                                              //
// Ce logiciel est un programme informatique faisant partie de la librairie     //
// SBCK (Statistical Bias Correction Kit). Cette librairie permet d'appliquer   //
// une correction de biais avec des méthodes non paramétriques, et propose      //
// diverses metrique entre Histograme Sparse en haute dimension.                //
//                                                                              //
// Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    //
// respectant les principes de diffusion des logiciels libres. Vous pouvez      //
// utiliser, modifier et/ou redistribuer ce programme sous les conditions       //
// de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     //
// sur le site "http://www.cecill.info".                                        //
//                                                                              //
// En contrepartie de l'accessibilité au code source et des droits de copie,    //
// de modification et de redistribution accordés par cette licence, il n'est    //
// offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    //
// seule une responsabilité restreinte pèse sur l'auteur du programme, le       //
// titulaire des droits patrimoniaux et les concédants successifs.              //
//                                                                              //
// A cet égard  l'attention de l'utilisateur est attirée sur les risques        //
// associés au chargement,  à l'utilisation,  à la modification et/ou au        //
// développement et à la reproduction du logiciel par l'utilisateur étant       //
// donné sa spécificité de logiciel libre, qui peut le rendre complexe à        //
// manipuler et qui le réserve donc à des développeurs et des professionnels    //
// avertis possédant  des  connaissances  informatiques approfondies.  Les      //
// utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       //
// logiciel à leurs besoins dans des conditions permettant d'assurer la         //
// sécurité de leurs systèmes et ou de leurs données et, plus généralement,     //
// à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           //
//                                                                              //
// Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       //
// pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    //
// termes.                                                                      //
//                                                                              //
//==============================================================================//
//==============================================================================//


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
