
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

//' @export
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
