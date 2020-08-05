
################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2019                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## This software is a computer program that is part of the SBCK (Statistical  ##
## Bias Correction Kit). This library makes it possible to perform bias       ##
## correction with non parametric methods, and give some metrics between      ##
## Sparse Histogram is high dimensions.                                       ##
##                                                                            ##
## This software is governed by the CeCILL-C license under French law and     ##
## abiding by the rules of distribution of free software.  You can  use,      ##
## modify and/ or redistribute the software under the terms of the CeCILL-C   ##
## license as circulated by CEA, CNRS and INRIA at the following URL          ##
## "http://www.cecill.info".                                                  ##
##                                                                            ##
## As a counterpart to the access to the source code and  rights to copy,     ##
## modify and redistribute granted by the license, users are provided only    ##
## with a limited warranty  and the software's author,  the holder of the     ##
## economic rights,  and the successive licensors  have only  limited         ##
## liability.                                                                 ##
##                                                                            ##
## In this respect, the user's attention is drawn to the risks associated     ##
## with loading,  using,  modifying and/or developing or reproducing the      ##
## software by the user in light of its specific status of free software,     ##
## that may mean  that it is complicated to manipulate,  and  that  also      ##
## therefore means  that it is reserved for developers  and  experienced      ##
## professionals having in-depth computer knowledge. Users are therefore      ##
## encouraged to load and test the software's suitability as regards their    ##
## requirements in conditions enabling the security of their systems and/or   ##
## data to be ensured and,  more generally, to use and operate it in the      ##
## same conditions as regards security.                                       ##
##                                                                            ##
## The fact that you are presently reading this means that you have had       ##
## knowledge of the CeCILL-C license and that you accept its terms.           ##
##                                                                            ##
################################################################################
################################################################################

################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2019                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## Ce logiciel est un programme informatique faisant partie de la librairie   ##
## SBCK (Statistical Bias Correction Kit). Cette librairie permet d'appliquer ##
## une correction de biais avec des méthodes non paramétriques, et propose    ##
## diverses metrique entre Histograme Sparse en haute dimension.              ##
##                                                                            ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et  ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez    ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions     ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA   ##
## sur le site "http://www.cecill.info".                                      ##
##                                                                            ##
## En contrepartie de l'accessibilité au code source et des droits de copie,  ##
## de modification et de redistribution accordés par cette licence, il n'est  ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le     ##
## titulaire des droits patrimoniaux et les concédants successifs.            ##
##                                                                            ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques      ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au      ##
## développement et à la reproduction du logiciel par l'utilisateur étant     ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à      ##
## manipuler et qui le réserve donc à des développeurs et des professionnels  ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les    ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du     ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la       ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,   ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         ##
##                                                                            ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez     ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les  ##
## termes.                                                                    ##
##                                                                            ##
################################################################################
################################################################################

# dataset_gaussian_exp_2d ##{{{

#' dataset_gaussian_exp_2d
#'
#' Generate a testing dataset such that the biased dataset is a distribution
#' of the the form Normal x Exp and the reference of the the form Exp x Normal.
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection) 
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_gaussian_exp_2d = function(n_samples)
{
	X0 = base::cbind( stats::rnorm(n_samples)             , stats::rexp( n_samples)  )
	Y0 = base::cbind( stats::rexp( n_samples)             , stats::rnorm(n_samples) )
	X1 = base::cbind( stats::rnorm(n_samples , mean = 5 ) , stats::rexp( n_samples)  )
	
	return( list( X0 = X0 , X1 = X1 , Y0 = Y0 ) )
}
##}}}

# dataset_gaussian_L_2d ##{{{

#' dataset_gaussian_L_2d
#'
#' Generate a testing dataset such that the biased dataset is a normal
#' distribution and reference a mixture a normal with a form in "L"
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection)
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_gaussian_L_2d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_gaussian_L_2d = function( n_samples )
{
	## Construction of X0 (biased period 0), X1 (biased period 1) and Y0 (reference period 0)
	size0  = as.integer(n_samples/2)
	size1  = n_samples - as.integer(n_samples/4)
	
	## Just a gaussian for X0
	X0 = ROOPSD::rmultivariate_normal( n = n_samples , mean = base::c(0.,0.) , cov = base::diag(2) )
	
	## A lightly complex gaussian for X1
	X1 = ROOPSD::rmultivariate_normal( n = n_samples , mean = base::c(1.,2.) , cov = matrix( base::c(2.,0,0,0.5) , nrow = 2 , ncol = 2 ) )
	
	## A very complex law for Y0
	Y0 = matrix( NA , nrow = n_samples , ncol = 2 )
	Y0[1:size0,]         = ROOPSD::rmultivariate_normal( n = size0             , mean = base::c(7.,7.)   , cov = matrix( base::c(2,0,0,0.5)   , nrow = 2 , ncol = 2 ) )
	Y0[size0:size1,]     = ROOPSD::rmultivariate_normal( n = size1 - size0 + 1 , mean = base::c(5.,9.)   , cov = matrix( base::c(0.5,0,0,2)   , nrow = 2 , ncol = 2 ) )
	Y0[size1:n_samples,] = ROOPSD::rmultivariate_normal( n = n_samples - size1 + 1 , mean = base::c(5.,12.5) , cov = matrix( base::c(0.2,0,0,0.2) , nrow = 2 , ncol = 2 ) )
	meanY0 = base::apply( Y0 , 2 , base::mean )
	meanX0 = base::apply( X0 , 2 , base::mean )
	diff = meanY0 - meanX0
	Y0 = base::t(base::apply( Y0 , 1 , function(x) { return(x - diff) } ))
	
	return( list( X0 = X0 , X1 = X1 , Y0 = Y0 ) )
}
##}}}

# dataset_gaussian_2d ##{{{

#' dataset_gaussian_2d
#'
#' Generate a testing dataset from random bivariate Gaussian distribution
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection)
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_gaussian_2d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_gaussian_2d = function(n_samples)
{
	CX0 = base::matrix( base::c(2.33,-0.38,-0.38,0.48) , nrow = 2 , ncol = 2 )
	X0  = ROOPSD::rmultivariate_normal( n_samples , mean = base::rep(0,2) , cov = CX0 )
	CX1 = base::matrix( base::c(0.1,-0.8,-0.8,3.37) , nrow = 2 , ncol = 2 )
	X1  = ROOPSD::rmultivariate_normal( n_samples , mean = base::rep(2,2) , cov = CX1 )
	CY0 = base::matrix( base::c(0.32,-0.18,-0.18,2.40) , nrow = 2 , ncol = 2 )
	Y0  = ROOPSD::rmultivariate_normal( n_samples , mean = base::rep(3,2) , cov = CY0 )
	
	return( list( X0 = X0 , X1 = X1 , Y0 = Y0 ) )
}
##}}}

# dataset_bimodal_reverse_2d ##{{{

#' dataset_bimodal_reverse_2d
#'
#' Generate a testing dataset from bimodale random bivariate Gaussian distribution
#'
#' @param n_samples [integer] numbers of samples drawn
#'        
#' @return [list] a list containing X0, X1 (biased in calibration/projection)
#'         and Y0 (reference in calibration)
#'
#' @examples
#' XY = SBCK::dataset_bimodal_reverse_2d(2000)
#' XY$X0 ## Biased in calibration period
#' XY$Y0 ## Reference in calibration period
#' XY$X1 ## Biased in projection period
#'
#' @export
dataset_bimodal_reverse_2d = function(n_samples)
{
	draw   = list( u = as.integer(n_samples/2) )
	draw[["l"]] = n_samples - draw$u
	lmY0   = list( u = base::c(5,-3)  , l = base::c(-3,3) )
	lcovY0 = list( u = 0.9 * diag(2)  , l = diag(2)       )
	lmX0   = list( u = base::c(0,0)   , l = base::c(2,2)  )
	lcovX0 = list( u = diag(2)        , l = 0.5 * diag(2) )
	lmX1   = list( u = base::c(-1,-1) , l = base::c(5,5)  )
	lcovX1 = list( u = 2 * diag(2)    , l = 0.1 * diag(2) )
	
	Y0 = NULL
	X0 = NULL
	X1 = NULL
	for( idx in base::c("u","l") )
	{
		Y0 = base::rbind( Y0 , ROOPSD::rmultivariate_normal( draw[[idx]] , lmY0[[idx]] , lcovY0[[idx]] ) )
		X0 = base::rbind( X0 , ROOPSD::rmultivariate_normal( draw[[idx]] , lmX0[[idx]] , lcovX0[[idx]] ) )
		X1 = base::rbind( X1 , ROOPSD::rmultivariate_normal( draw[[idx]] , lmX1[[idx]] , lcovX1[[idx]] ) )
	}
	
	return( list( X0 = X0 , X1 = X1 , Y0 = Y0 ) )
}

##}}}

