
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

#' Univariate Random Variable from histogram
#'
#' Build an univariate random variable from an histogram of a dataset X.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param bins [vector of NULL]
#'        A vector of bins.
#'        If NULL, it is estimating.
#' @param X [vector]
#'        Vector of data
#' @param q [vector]
#'        Vector of quantiles
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(X,bins)}}{This method is used to create object of this class with \code{rv_histogram}}
#'   \item{\code{rvs(size)}}{Random values generator from histogram estimated}.
#'   \item{\code{cdf(X)}}{Cumulative Distribution Function.}.
#'   \item{\code{icdf(X)}}{Inverse of Cumulative Distribution Function.}.
#'   \item{\code{sf(X)}}{Survival Function (1-CDF).}.
#'   \item{\code{isf(X)}}{Inverse of Survival Function.}.
#' }
#' @examples
#' ## Realizations of a random variable
#' X = stats::rnorm( 10000 )
#'
#' ## Estimation of random variable
#' rvX = SBCK::rv_histogram$new(X)
#' 
#' ## cdf and sf
#' x = base::seq( -2 , 2 , 0.001 )
#' cdfx = rvX$cdf(x)
#' sfx = rvX$sf(x)
#'
#' ## icdf and isf
#' q = base::seq( 0 , 1 , 0.001 )
#' icdfq = rvX$icdf(q)
#' isfq = rvX$isf(q)
#'
#' @export
rv_histogram = R6::R6Class( "rv_histogram" ,
	public = list(
	
	###############
	## Arguments ##
	###############
	
	bins  = NULL,
	min   = NULL,
	max   = NULL,
	p     = NULL,
	c     = NULL,
	nbins = NULL,
	
	#################
	## Constructor ##
	#################
	
	initialize = function( X , bins = NULL )
	{
		self$min = base::min(X)
		self$max = base::max(X)
		
		## Bins construction
		bin_width = 0
		if( is.null(bins) )
		{
			bin_width = bin_width_estimator(X)
			self$bins = base::seq( self$min - bin_width , self$max + bin_width , bin_width )
		}
		else
		{
			bin_width = min( diff(bins) )
			self$bins = bins
		}
		
		## Histogram
		hist = graphics::hist( X , breaks = self$bins , plot = FALSE )
		self$p = hist$density / base::sum(hist$density)
		self$c = hist$mids
		self$nbins = length(self$p)
		
		## CDF and iCDF function
		private$cdffn = stats::ecdf(X)
		x = base::seq( self$min - bin_width , self$max + bin_width , 1e-3 * bin_width )
		quants = private$cdffn(x)
		private$icdffn = stats::approxfun( quants , x , yleft = self$min - bin_width , yright = self$max + bin_width )
	},
	
	rvs = function( size )
	{
		idx = sample( 1:self$nbins , size = size , replace = TRUE , prob = self$p )
		invisible(self$c[idx])
	},
	
	cdf = function( X )
	{
		invisible(private$cdffn(X))
	},
	
	icdf = function( q )
	{
		invisible(private$icdffn(q))
	},
	
	sf = function( X )
	{
		invisible( 1. - private$cdffn(X) )
	},
	
	isf = function( q )
	{
		invisible(private$icdffn(1. - q))
	}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	cdffn = NULL,
	icdffn = NULL
	)
)
