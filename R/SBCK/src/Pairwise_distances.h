
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


#ifndef SBCK_PAIRWISE_DISTANCES_INCLUDED
#define SBCK_PAIRWISE_DISTANCES_INCLUDED

#include <iostream>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <cmath>

//---------//
// metrics //
//---------//

typedef Rcpp::NumericMatrix::Row    MatRow ;
typedef Rcpp::NumericMatrix::Column MatCol ;
typedef std::function<double(MatRow,MatRow)> MetricType ;

inline double sqeuclidean_metric( MatRow x , MatRow y ) //{{{
{
	double dist(0) ;
	std::size_t size(x.size()) ;
	for( std::size_t s = 0 ; s < size ; ++s )
		dist += std::pow( x[s] - y[s] , 2 ) ;
	return dist ;
} //}}}

inline double euclidean_metric( MatRow x , MatRow y ) //{{{
{
	return std::sqrt( sqeuclidean_metric( x , y ) ) ;
} //}}}

inline double chebyshev_metric( MatRow x , MatRow y ) //{{{
{
	double dist(-1) ;
	std::size_t size(x.size()) ;
	for( std::size_t s = 0 ; s < size ; ++s )
		dist = std::max( dist , std::abs( x[s] - y[s] ) ) ;
	return dist ;
} //}}}

inline double logeuclidean_metric( MatRow x , MatRow y ) //{{{
{
	double dist(0) ;
	std::size_t size(x.size()) ;
	for( std::size_t s = 0 ; s < size ; ++s )
		dist += std::pow( x[s] - y[s] , 2 ) ;
	return std::log( dist ) / 2. ;
} //}}}

inline MetricType metric_choosen( std::string str_metric ) //{{{
{
	if( str_metric == "sqeuclidean" )
		return sqeuclidean_metric ;
	else if( str_metric == "chebyshev" )
		return chebyshev_metric ;
	else if( str_metric == "logeuclidean" )
		return logeuclidean_metric ;
	else
		return euclidean_metric ;
} //}}}


//----------//
// Pairwise //
//----------//

// cpp_pairwise_distances_XYstr {{{

//' cpp_pairwise_distances_XYstr
//' 
//' Pairwise distances between two differents matrix X and Y with a 
//' compiled str_metric. DO NOT USE, use CDSK::pairwise_distances
//'
//' @usage cpp_pairwise_distances_XYstr(X,Y,str_metric)
//' @param X [Rcpp::NumericMatrix] Matrix
//' @param Y [Rcpp::NumericMatrix] Matrix
//' @param str_metric [std::string] c++ string
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericMatrix cpp_pairwise_distances_XYstr( Rcpp::NumericMatrix X , Rcpp::NumericMatrix Y , std::string str_metric )
{
	// Quelques paramètres
	std::size_t size0 = X.nrow() ;
	std::size_t size1 = Y.nrow() ;
	MetricType metric = metric_choosen( str_metric ) ;
	
	// Calcul des distances
	Rcpp::NumericMatrix Dist( size0 , size1 ) ;
	for( std::size_t i = 0 ; i < size0 ; ++i )
	{
		for( std::size_t j = 0 ; j < size1 ; ++j )
		{
			Dist(i,j) = metric( X(i,Rcpp::_) , Y(j,Rcpp::_) ) ;
		}
	}
	return Dist ;
} //}}}

// cpp_pairwise_distances_Xstr {{{

//' cpp_pairwise_distances_Xstr
//' 
//' Pairwise distances between X and themselves with a compiled 
//' str_metric. DO NOT USE, use CDSK::pairwise_distances
//'
//' @usage cpp_pairwise_distances_Xstr(X,str_metric)
//' @param X [Rcpp::NumericMatrix] Matrix
//' @param str_metric [std::string] c++ string
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericMatrix cpp_pairwise_distances_Xstr( Rcpp::NumericMatrix X , std::string str_metric )
{
	// Quelques paramètres
	std::size_t size = X.nrow() ;
	MetricType metric = metric_choosen( str_metric ) ;
	
	// Calcul des distances
	Rcpp::NumericMatrix Dist( size , size ) ;
	for( std::size_t i = 0 ; i < size ; ++i )
	{
		for( std::size_t j = i ; j < size ; ++j )
		{
			Dist(i,j) = metric( X(i,Rcpp::_) , X(j,Rcpp::_) ) ;
			Dist(j,i) = Dist(i,j) ;
		}
	}
	return Dist ;
} //}}}

// cpp_pairwise_distances_XYCall {{{

//' cpp_pairwise_distances_XYCall
//' 
//' Pairwise distances between X  and Y with a R function (metric).
//' DO NOT USE, use CDSK::pairwise_distances
//'
//' @usage cpp_pairwise_distances_XYCall(X,Y,metric)
//' @param X [Rcpp::NumericMatrix] Matrix
//' @param Y [Rcpp::NumericMatrix] Matrix
//' @param metric [Rcpp::Function] R function
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericMatrix cpp_pairwise_distances_XYCall( Rcpp::NumericMatrix X , Rcpp::NumericMatrix Y , Rcpp::Function metric )
{
	// Quelques paramètres
	std::size_t size0 = X.nrow() ;
	std::size_t size1 = Y.nrow() ;
	
	// Calcul des distances
	Rcpp::NumericMatrix Dist( size0 , size1 ) ;
	for( std::size_t i = 0 ; i < size0 ; ++i )
	{
		for( std::size_t j = 0 ; j < size1 ; ++j )
		{
			Dist(i,j) = Rcpp::as<double>( metric( X(i,Rcpp::_) , Y(j,Rcpp::_) ) ) ;
		}
	}
	return Dist ;
} //}}}

// cpp_pairwise_distances_XCall {{{

//' cpp_pairwise_distances_XCall
//' 
//' Pairwise distances between X and themselves with a R function (metric).
//' DO NOT USE, use CDSK::pairwise_distances
//'
//' @usage cpp_pairwise_distances_XCall(X,metric)
//' @param X [Rcpp::NumericMatrix] Matrix
//' @param metric [Rcpp::Function] R function
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericMatrix cpp_pairwise_distances_XCall( Rcpp::NumericMatrix X , Rcpp::Function metric )
{
	// Quelques paramètres
	std::size_t size = X.nrow() ;
	
	// Calcul des distances
	Rcpp::NumericMatrix Dist( size , size ) ;
	for( std::size_t i = 0 ; i < size ; ++i )
	{
		for( std::size_t j = i ; j < size ; ++j )
		{
			Dist(i,j) = Rcpp::as<double>(metric( X(i,Rcpp::_) , X(j,Rcpp::_) )) ;
			Dist(j,i) = Dist(i,j) ;
		}
	}
	return Dist ;
} //}}}

#endif
