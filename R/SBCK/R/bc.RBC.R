
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

#' RBC (Random Bias Correction) method
#'
#' @description
#' Perform a multivariate bias correction of X with respect to Y randomly.
#'
#' @details
#' Only for comparison.
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#'
#' ## Bias correction
#' ## Step 1 : construction of the class RBC
#' rbc = SBCK::RBC$new() 
#' ## Step 2 : Fit the bias correction model
#' rbc$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = rbc$predict(X1,X0) 
#' ## Z$Z0 # BC of X0
#' ## Z$Z1 # BC of X1
#' @export
RBC = R6::R6Class( "RBC" ,
	
	
	public = list(
	
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new RBC object.
	#'
	#' @return A new `RBC` object.
	initialize = function()
	{},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection, can be
    #'        NULL for stationary BC method
    #' @return NULL
	fit = function( Y0 , X0 , X1 = NULL )
	{
		private$Y0 = if( !is.matrix(Y0) ) matrix( Y0 , nrow = length(Y0) , ncol = 1 ) else Y0
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction. Use named keywords to use stationary or
    #' non-stationary method.
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features or NULL] Model in projection
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 = NULL , X0 = NULL )
	{
		Z0 = NULL
		Z1 = NULL
		if( !is.null(X0) )
		{
			if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
			idx = base::sample( 1:base::nrow(private$Y0) , base::nrow(X0) , replace = TRUE )
			Z0 = private$Y0[idx,]
		}
		
		if( !is.null(X1) )
		{
			if( !is.matrix(X1) ) X1 = matrix( X1 , nrow = length(X0) , ncol = 1 )
			idx = base::sample( 1:base::nrow(private$Y0) , base::nrow(X1) , replace = TRUE )
			Z1 = private$Y0[idx,]
		}
		
		if( !is.null(X0) && !is.null(X1) )
		{
			return( list( Z0 = Z0 , Z1 = Z1 ) )
		}
		else if( is.null(X0) )
		{
			return(Z1)
		}
		else
		{
			return(Z0)
		}
	}
	##}}}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	Y0 = NULL
	
	)
)
