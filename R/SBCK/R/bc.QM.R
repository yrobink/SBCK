
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

#' Quantile Mapping method
#'
#' Perform an univariate bias correction of X with respect to Y (correction is applied margins by margins).
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param bins [list of vector of NULL]
#'        A list of bins for each margins.
#'        If NULL, it is estimating during the fit
#' @param Y  [matrix]
#'        A matrix containing references (nrow = n_samples, ncol = n_features)
#' @param X [matrix]
#'        A matrix containing biased data (nrow = n_samples, ncol = n_features)
#'
#' @return Object of \code{\link{R6Class}} with methods for bias correction
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(bins)}}{This method is used to create object of this class with \code{QM}}
#'   \item{\code{fit(Y,X)}}{Fit the bias correction model from Y and X}.
#'   \item{\code{predict(X)}}{Perform the bias correction of X with respect to Y.}.
#' }
#' @references Panofsky, H. A. and Brier, G. W.: Some applications of statistics to meteorology, Mineral Industries Extension Services, College of Mineral Industries, Pennsylvania State University, 103 pp., 1958.
#' @references Wood, A. W., Leung, L. R., Sridhar, V., and Lettenmaier, D. P.: Hydrologic Implications of Dynamical and Statistical Approaches to Downscaling Climate Model Outputs, Clim. Change, 62, 189–216, https://doi.org/10.1023/B:CLIM.0000013685.99609.9e, 2004.
#' @references Déqué, M.: Frequency of precipitation and temperature extremes over France in an anthropogenic scenario: Model results and statistical correction according to observed values, Global Planet. Change, 57, 16–26, https://doi.org/10.1016/j.gloplacha.2006.11.030, 2007.
#' @examples
#' ## Two bivariate random variables (rnorm and rexp are inverted between ref and bias)
#' Y = base::cbind( rnorm(10000) , rexp(10000)  )
#' X = base::cbind( rexp(10000)  , rnorm(10000) )
#'
#' ## Bias correction
#' ## Step 1 : construction of the class QM 
#' qm = SBCK::QM$new() 
#' ## Step 2 : Fit the bias correction model
#' qm$fit( Y , X )
#' ## Step 3 : perform the bias correction, uX is the correction of
#' ## X with respect to the estimation of Y
#' uX = qm$predict(X) 
#'
#' @export
QM = R6::R6Class( "QM" ,
	
	public = list(
	
	
	###############
	## Arguments ##
	###############
	
	n_features = 0,
	bins = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( bins = NULL )
	{
		self$bins = bins
	},
	
	fit = function( Y , X )
	{
		## Dimension and data formating
		if( class(Y) == "numeric" )
			Y = matrix( Y , nrow = length(Y) , ncol = 1 )
		if( class(X) == "numeric" )
			X = matrix( X , nrow = length(X) , ncol = 1 )
		
		self$n_features = dim(X)[2]
		
		## Bins
		if( is.null(self$bins) )
		{
			bin_width = SBCK::common_bin_width_estimator( list(X,Y) )
			self$bins = list()
			for( i in 1:self$n_features )
			{
				self$bins[[i]] = base::seq( min(X[,i],Y[,i]) - bin_width[i] , max(X[,i],Y[,i]) + bin_width[i] , bin_width[i] )
			}
		}
		
		## Random variable
		for( i in 1:self$n_features )
		{
			private$rvX[[i]] = rv_histogram$new( X[,i] , self$bins[[i]] )
			private$rvY[[i]] = rv_histogram$new( Y[,i] , self$bins[[i]] )
		}
	},
	
	predict = function( X )
	{
		if( class(X) == "numeric" )
		{
			X = matrix( X , nrow = length(X) , ncol = 1 )
		}
		
		Z = matrix( NA , nrow = dim(X)[1] , ncol = self$n_features )
		
		for( i in 1:self$n_features )
		{
			Z[,i] = private$rvY[[i]]$icdf( private$rvX[[i]]$cdf( X[,i] ) )
		}
		return(Z)
	}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	rvX = list(),
	rvY = list()
	)
)
