
##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the SBCK (Statistical    ##
## Bias Correction Kit). This library makes it possible to perform bias         ##
## correction with non parametric methods, and give some metrics between Sparse ##
## Histogram is high dimensions.                                                ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## SBCK (Statistical Bias Correction Kit). Cette librairie permet d'appliquer   ##
## une correction de biais avec des méthodes non paramétriques, et propose      ##
## diverses metrique entre Histograme Sparse en haute dimension.                ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################

#' bin_width_estimator method
#'
#' Lenght of cell to compute an histogram
#'
#' @param X [matrix]
#'        A matrix containing data, nrow = n_samples, ncol = n_features
#' @param method [string]
#'        Method to estimate bin_width, values are "auto", "FD" (Friedman Draconis,
#'        robust over outliners) or "Sturges". If "auto" is used and if nrow(X) < 1000,
#'        "Sturges" is used, else "FD" is used.
#'        
#' @return [vector] Lenght of bins
#'
#' @importFrom stats quantile
#' @examples
#' X = base::cbind( stats::rnorm( n = 2000 ) , stats::rexp(2000) )
#' binw_width = SBCK::bin_width_estimator( X , method = "auto" ) ## Friedman Draconis is used
#' X = stats::rnorm( n = 500 )
#' binw_width = SBCK::bin_width_estimator( X , method = "auto" ) ## Sturges is used
#' 
#' @export
bin_width_estimator = function( X , method = "auto" ) 
{
	if( class(X) == "numeric" )
		X = matrix( X , nrow = length(X) , ncol = 1 )
	
	n_samples  = dim(X)[1]
	n_features = dim(X)[2]
	
	## Find method to use
	if( method == "auto" && n_samples < 1000 )
	{
		method = "Sturges"
	}
	
	## Find bin_width
	bin_width = rep( 0 , n_features )
	if( method == "Sturges" )
	{
		nh = log2(n_samples) + 1
		bin_width = rep( 1. / nh , n_features )
	}
	else ## FD (Freedman Diaconis) method, robust over outliers
	{
		pow = n_samples^(1./3.)
		q    = base::apply( X , 2 , quantile , probs = base::c( 0.25 , 0.75 ) )
		bin_width = 2 * ( q[2,] - q[1,] ) / pow
	}
	
	invisible(as.vector(bin_width))
}


#' common_bin_width_estimator method
#'
#' Common lenght of cell to compute an histogram from many dataset
#'
#' @param lX [list of matrix]
#'        A list of matrix containing data, nrow = n_samples, ncol = n_features
#' @param method [string]
#'        Method to estimate bin_width, values are "auto", "FD" (Friedman Draconis,
#'        robust over outliners) or "Sturges". If "auto" is used and if nrow(X) < 1000,
#'        "Sturges" is used, else "FD" is used.
#'        
#' @return [vector] Lenght of bins
#'
#' @examples
#' X = base::cbind( stats::rnorm( n = 2000 ) , stats::rexp(2000) )
#' Y = base::cbind( stats::rnorm( n = 2000 ) , stats::rexp(2000) )
#' Z = base::cbind( stats::rnorm( n = 2000 ) , stats::rexp(2000) )
#' bin_width = common_bin_width_estimator( list(X,Y,Z) )
#'
#' @export
common_bin_width_estimator = function( lX , method = "auto" )
{
	bw = matrix( NA , nrow = length(lX) , ncol = dim(lX[[1]])[2] )
	for( i in 1:length(lX) )
		bw[i,] = SBCK::bin_width_estimator( lX[[i]] , method )
	invisible( base::apply( bw , 2 , base::min ) )
}

