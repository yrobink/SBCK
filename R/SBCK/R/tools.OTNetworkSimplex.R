
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


#' Optimal Transport Network Simplex solver
#'
#' @description
#' Solve the optimal transport problem 
#'
#' @details
#' use the network simplex algorithm
#'
#' @references Bazaraa, M. S., Jarvis, J. J., and Sherali, H. D.: Linear
#'             Programming and Network Flows, 4th edn., John Wiley & Sons, 2009.
#'
#' @examples
#' ## Define two dataset
#' X = stats::rnorm(2000)
#' Y = stats::rnorm(2000 , mean = 5 )
#' bw = base::c(0.1)
#' muX = SBCK::SparseHist( X , bw )
#' muY = SBCK::SparseHist( Y , bw )
#' 
#' ## Find solution
#' ot = OTNetworkSimplex$new()
#' ot$fit( muX , muY )
#' 
#' print( sum(ot$plan) ) ## Must be equal to 1
#' print( ot$success )   ## If solve is success
#' print( sqrt(sum(ot$plan * ot$C)) ) ## Cost of plan
#'
#' @export
OTNetworkSimplex = R6::R6Class( "OTNetworkSimplex" ,
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field p [double] Power of the plan
	p       = NULL,
	#' @field plan [matrix] transport plan
	plan    = NULL,
	#' @field success [bool] If the fit is a success or not
	success = NULL,
	#' @field C [matrix] Cost matrix
	C       = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description 
    #' Create a new OTNetworkSimplex object.
	#' @param p [double] Power of the plan
	#'
	#' @return A new `OTNetworkSimplex` object.
	initialize = function( p = 2 )
	{
		self$p = p
	},
	##}}}
	
	## fit ##{{{
	#' @description
	#' Fit the OT plan
	#' @param muX0 [SparseHist or OTHist] Source histogram to move
	#' @param muX1 [SparseHist or OTHist] Target histogram
	#' @param C [matrix or NULL] Cost matrix (without power p) between muX0 and
	#'       muX1, if NULL pairwise_distances is called with Euclidean distance.
	#'
	#' @return NULL
	fit = function( muX0 , muX1 , C = NULL )
	{
		self$C = C
		if( is.null(self$C) )
			self$C = pmetric::pairwise_distances( muX0$c , muX1$c )^self$p
		out = SBCK::network_simplex( muX0$p , muX1$p , self$C )
		self$plan    = out$plan
		self$success = out$success == TRUE
	}
	##}}}
	
	)
	
	
	######################
	## Private elements ##
	######################
	
)
