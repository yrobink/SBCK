
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
//' compiled str_metric. DO NOT USE, use SBCK::pairwise_distances
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
//' str_metric. DO NOT USE, use SBCK::pairwise_distances
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
//' DO NOT USE, use SBCK::pairwise_distances
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
//' DO NOT USE, use SBCK::pairwise_distances
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
