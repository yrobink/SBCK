
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


#' Slope stopping criteria
#'
#' @description
#' Class which send a stop signal when a time series stay constant.
#'
#' @details
#' Test the slope.
#'
#' @examples
#' stop_slope = SlopeStoppingCriteria$new( 20 , 500 , 1e-3 )
#' x = 0
#' while(!stop_slope$stop)
#' {
#' 	stop_slope$append(base::exp(-x))
#' 	x = x + 0.1
#' }
#' print(stop_slope$nit)
#' @export
SlopeStoppingCriteria = R6::R6Class( "SlopeStoppingCriteria" ,
	
	public = list(
	
	################
	## Parameters ##
	################
	
	#' @field minit [integer] Minimal number of iterations. At least 3.
	minit    = NULL,
	#' @field maxit [integer] Maximal number of iterations.
	maxit    = NULL,
	#' @field nit [integer] Number of iterations.
	nit      = NULL,
	#' @field tol [float] Tolerance to control if slope is close to zero
	tol      = NULL,
	#' @field stop [bool] If we stop
	stop     = NULL,
	#' @field criteria [vector] State of criteria
	criteria = NULL,
	#' @field slope [vector] Values of slope
	slope    = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new SlopeStoppingCriteria object.
	#' @param minit [integer] Minimal number of iterations. At least 3.
	#' @param maxit [integer] Maximal number of iterations.
	#' @param tol [float] Tolerance to control if slope is close to zero
	#'
	#' @return A new `SlopeStoppingCriteria` object.
	initialize = function( minit , maxit , tol )
	{
		self$minit = max(minit,3)
		self$maxit = maxit
		self$tol = tol
		self$nit = 1
		self$stop = FALSE
		self$criteria = numeric(maxit)
		self$slope = numeric(maxit)
	},
	##}}}
	
	## reset ##{{{
	#' @description
    #' Reset the class
	#'
	#' @return NULL
	reset = function()
	{
		self$nit = 1
		self$stop = FALSE
		self$criteria = numeric(self$maxit)
		self$slope = numeric(self$maxit)
	},
	##}}}
	
	## append ##{{{
	#' @description
    #' Add a new value
    #' @param value [double] New metrics
	#'
	#' @return NULL
	append = function(value)
	{
		self$criteria[self$nit] = value
		if( self$nit >= 2 )
		{
			x = 1:self$nit
			y = self$criteria[1:self$nit]
			lm = stats::lm( y ~ x )
			self$slope[self$nit-1] = as.vector(lm$coefficients[2])
		}
		self$stop = (self$nit >= self$minit) && (self$nit >= self$maxit - 1 || base::abs(self$slope[self$nit-1]) < self$tol)
		self$nit = self$nit + 1
	}
	##}}}
	
	)
)


