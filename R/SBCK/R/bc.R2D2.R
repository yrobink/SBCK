
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

#' R2D2 (Rank Resampling for Distributions and Dependences) method
#'
#' @description
#' Perform a multivariate (non stationary) bias correction.
#'
#' @details
#' Use rankshuffle in calibration and projection period with CDFt
#'
#' @references Vrac, M.: Multivariate bias adjustment of high-dimensional
#'             climate simulations: the Rank Resampling for Distributions and
#'             Dependences (R2 D2 ) bias correction, Hydrol. Earth Syst. Sci.,
#'             22, 3175–3196, https://doi.org/10.5194/hess-22-3175-2018, 2018.
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#' ## Bias correction
#' ## Step 1 : construction of the class R2D2 
#' r2d2 = SBCK::R2D2$new() 
#' ## Step 2 : Fit the bias correction model
#' r2d2$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction
#' Z = r2d2$predict(X1,X0) 
#'
#' @export
R2D2 = R6::R6Class( "R2D2" ,
	
	inherit = SBCK::CDFt,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field irefs [vector of int] Indexes for shuffle. Defaults is base::c(1)
	irefs = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new R2D2 object.
	#' @param irefs [vector of int] Indexes for shuffle. Defaults is base::c(1)
	#'        model
	#' @param ... [] all others arguments are passed to CDFt class.
	#'
	#' @return A new `R2D2` object.
	initialize = function( irefs = base::c(1) , ... )
	{
		kwargs = list(...)
		base::do.call( super$initialize , kwargs )
		self$irefs = irefs
		private$ssr = SBCK::SchaakeShuffleRef$new( 1 )
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return NULL
	fit = function( Y0 , X0 , X1 )
	{
		super$fit( Y0 , X0 , X1 )
		private$ssr$fit(Y0)
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 , X0 = NULL )
	{
		Z = super$predict(X1,X0)
		if( is.null(X0) )
		{
			Z1 = base::array( NA , dim = base::c( base::dim(Z) , length(self$irefs) ) )
			for( i in 1:length(self$irefs) )
			{
				private$ssr$ref = self$irefs[i]
				Z1[,,i] = private$ssr$predict(Z)
			}
			if( length(self$irefs) == 1 )
				return(Z1[,,1])
			return(Z1)
		}
		
		Z0 = base::array( NA , dim = base::c( base::dim(Z$Z0) , length(self$irefs) ) )
		Z1 = base::array( NA , dim = base::c( base::dim(Z$Z1) , length(self$irefs) ) )
		for( i in 1:length(self$irefs) )
		{
			private$ssr$ref = self$irefs[i]
			Z0[,,i] = private$ssr$predict(Z$Z0)
			Z1[,,i] = private$ssr$predict(Z$Z1)
		}
		if( length(self$irefs) == 1 )
		{
			Z0 = Z0[,,1]
			Z1 = Z1[,,1]
		}
		
		return( list( Z1 = Z1 , Z0 = Z0 ) )
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
	
	ssr  = NULL
	)
)
