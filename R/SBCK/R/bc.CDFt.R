
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

##################################################################################
##################################################################################
##                                                                              ##
## Original author  : Mathieu Vrac                                              ##
## Contact          : mathieu.vrac@lsce.ipsl.fr                                 ##
##                                                                              ##
## Notes   : CDFt is the re-implementation of the function CDFt of R package    ##
##           "CDFt" developped by Mathieu Vrac, available at                    ##
##           https://cran.r-project.org/web/packages/CDFt/index.html            ##
##           This code is governed by the CeCILL-C license with the             ##
##           authorization of Mathieu Vrac                                      ##
##                                                                              ##
##################################################################################
##################################################################################

###############
## Libraries ##
###############

###############
## Functions ##
###############


## CDFt Correction method

#' CDFt method (Cumulative Distribution Function transfer)
#'
#' Perform an univariate bias correction of X with respect to Y (correction is applied margins by margins).
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param ...
#'        Many named arguments listed below
#' @param distX0 [A ROOPSD:: distribution or a list of them]
#'        Describe the law of each margins. A list permit to use different laws for each margins. Default is rv_histogram.
#' @param distY0 [A ROOPSD:: distribution or a list of them]
#'        Describe the law of each margins. A list permit to use different laws for each margins. Default is rv_histogram.
#' @param distX1 [A ROOPSD:: distribution or a list of them]
#'        Describe the law of each margins. A list permit to use different laws for each margins. Default is rv_histogram.
#' @param distY1 [A ROOPSD:: distribution or a list of them]
#'        Describe the law of each margins. A list permit to use different laws for each margins. Default is rv_histogram.
#' @param n_features  [NULL or integer]
#'        Normaly infered during fit, but if distX0, distX1 and distY0 are simultaneously frozen, must be set during initialization.
#' @param Y0  [matrix]
#'        A matrix containing references during calibration period
#' @param X0 [matrix]
#'        A matrix containing biased data during calibration period
#' @param X1 [matrix]
#'        A matrix containing biased data during projection period
#'
#' @return Object of \code{\link{R6Class}} with methods for bias correction
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(...)}}{This method is used to create object of this class with \code{CDFt}}
#'   \item{\code{fit(Y0,X0,X1)}}{Fit the bias correction model from Y0 and X0 and X1}.
#'   \item{\code{predict(X1,X0)}}{Perform the bias correction of X1 with respect to the estimaton of Y1.}.
#' }
#' @references Michelangeli, P.-A., Vrac, M., and Loukos, H.: Probabilistic downscaling approaches: Application to wind cumulative distribution functions, Geophys. Res. Lett., 36, L11708, https://doi.org/10.1029/2009GL038401, 2009.
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#' ## Bias correction
#' ## Step 1 : construction of the class CDFt 
#' cdft = SBCK::CDFt$new() 
#' ## Step 2 : Fit the bias correction model
#' cdft$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction, Z is a list containing
#' ## corrections
#' Z = cdft$predict(X1,X0) 
#' Z$Z0 ## Correction in calibration period
#' Z$Z1 ## Correction in projection period
#' @export
CDFt = R6::R6Class( "CDFt" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	n_features = 0,
	tol        = 1e-3,
	
	distY0 = NULL,
	distY1 = NULL,
	distX0 = NULL,
	distX1 = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function(...) ##{{{
	{
		kwargs = list(...)
		self$distY0 = DistHelper$new( dist = kwargs[["distY0"]] , kwargs = kwargs[["kwargsY0"]] )
		self$distY1 = DistHelper$new( dist = kwargs[["distY1"]] , kwargs = kwargs[["kwargsY1"]] )
		self$distX0 = DistHelper$new( dist = kwargs[["distX0"]] , kwargs = kwargs[["kwargsX0"]] )
		self$distX1 = DistHelper$new( dist = kwargs[["distX1"]] , kwargs = kwargs[["kwargsX1"]] )
		self$n_features = kwargs[["n_features"]]
		self$tol        = if( !is.null(kwargs[["tol"]]) ) kwargs[["tol"]] else 1e-3
	},
	##}}}
	
	fit = function( Y0 , X0 , X1 )##{{{
	{
		## Dimension and data formating
		##=============================
		if( class(Y0) == "numeric" ) Y0 = matrix( Y0 , nrow = length(Y0) , ncol = 1 )
		if( class(X0) == "numeric" ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		if( class(X1) == "numeric" ) X1 = matrix( X1 , nrow = length(X1) , ncol = 1 )
		
		## Set features
		##=============
		if( is.null(self$n_features ) )
		{
			if( !is.null(Y0) )
				self$n_features = base::ncol(Y0)
			else if( !is.null(X0) )
				self$n_features = base::ncol(X0)
			else if( !is.null(X1) )
				self$n_features = base::ncol(X1)
			else
				base::stop("SBCK::CDFt: If Y0 == X0 == X1 == NULL, n_features must be set during initialization")
		}
		self$distY0$set_features(self$n_features)
		self$distY1$set_features(self$n_features)
		self$distX0$set_features(self$n_features)
		self$distX1$set_features(self$n_features)
		
		
		## Start fit itself
		##=================
		for( i in 1:self$n_features )
		{
			self$distY0$fit( Y0[,i] , i )
			self$distX0$fit( X0[,i] , i )
			self$distX1$fit( X1[,i] , i )
			
			if( self$distY1$is_frozen(i) )
			{
				self$distY1$law[[i]] = self$distY1$dist[[i]]
			}
			else
			{
				Y1 = NULL
				if( self$distY0$is_parametric(i) && self$distX0$is_parametric(i) && self$distX1$is_parametric(i) )
				{
					Y1 = self$distX1$law[[i]]$icdf( self$distX0$law[[i]]$cdf( self$distY0$law[[i]]$icdf( self$distX1$law[[i]]$cdf( X1[,i] ) ) ) )
				}
				else
				{
					Y0uni = if( is.null(Y0) ) self$distY0$law[[i]]$rvs(10000) else Y0[,i]
					X0uni = if( is.null(X0) ) self$distX0$law[[i]]$rvs(10000) else X0[,i]
					X1uni = if( is.null(X1) ) self$distX1$law[[i]]$rvs(10000) else X1[,i]
					Y1 = private$infer_Y1( Y0uni , X0uni , X1uni , i ) ## infer Y1 here
				}
				self$distY1$fit( Y1 , i )
			}
		}
	},
	##}}}
	
	predict = function( X1 , X0 = NULL )##{{{
	{
		if( class(X0) == "numeric" ) X0 = matrix( X0 , nrow = length(X1) , ncol = 1 )
		if( class(X1) == "numeric" ) X1 = matrix( X1 , nrow = length(X1) , ncol = 1 )
		
		if( is.null(X0) )
		{
			Z1 = matrix( NA , nrow = base::nrow(X1) , ncol = base::ncol(X1) )
			for( i in 1:self$n_features )
			{
				cdf = self$distX1$law[[i]]$cdf( X1[,i] )
				cdf[!(cdf < 1)] = 1-self$tol
				cdf[!(cdf > 0)] = self$tol
				Z1[,i] = self$distY1$law[[i]]$icdf( cdf )
			}
			return(Z1)
		}
		else
		{
			Z0 = matrix( NA , nrow = base::nrow(X0) , ncol = base::ncol(X0) )
			Z1 = matrix( NA , nrow = base::nrow(X1) , ncol = base::ncol(X1) )
			for( i in 1:self$n_features )
			{
				cdf = self$distX0$law[[i]]$cdf( X0[,i] )
				cdf[!(cdf < 1)] = 1-self$tol
				cdf[!(cdf > 0)] = self$tol
				Z0[,i] = self$distY0$law[[i]]$icdf( cdf )
				
				cdf = self$distX1$law[[i]]$cdf( X1[,i] )
				cdf[!(cdf < 1)] = 1-self$tol
				cdf[!(cdf > 0)] = self$tol
				Z1[,i] = self$distY1$law[[i]]$icdf( cdf )
			}
			return( list( Z1 = Z1 , Z0 = Z0 ) )
		}
		
		return(Z1)
	}##}}}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(##{{{
	
	###############
	## Arguments ##
	###############
	
	qmX1Y1 = list(),
	qmX0Y0 = NULL,
	diff = NULL,
	
	
	#############
	## Methods ##
	#############
	
	infer_Y1 = function( Y0 , X0 , X1 , idx )##{{{
	{
		mY0 = base::mean(Y0)
		mX0 = base::mean(X0)
		
		X0s = X0 + mY0 - mX0
		X1s = X1 + mY0 - mX0
		
		rvY0 = self$distY0$law[[idx]]
		rvX0 = base::do.call( self$distX0$dist[[idx]]$new , self$distX0$kwargs )
		rvX0$fit(X0s)
		rvX1 = base::do.call( self$distX1$dist[[idx]]$new , self$distX1$kwargs )
		rvX1$fit(X1s)
		
		xdiff = base::abs( base::mean(X1) - base::mean(X0) )
		xdev = 2
		dev_ok = FALSE
		
		while( !dev_ok )
		{
			dev_ok = TRUE
			xmin = base::min(Y0,X0,X1) - xdev * xdiff
			xmax = base::max(Y0,X0,X1) + xdev * xdiff
			
			x = base::seq( xmin , xmax , length = 200 )
			cdfY0 = rvY0$cdf(x)
			cdfX0 = rvX0$cdf(x)
			cdfX1 = rvX1$cdf(x)
			cdfY1 = rvY0$cdf(rvX0$icdf(cdfX1))
			
			if( base::min(Y0) < base::min(X1s) )
			{
				i = base::max( base::which( x < rvY0$icdf(cdfY1[1]) ) )
				j = base::max( base::which( x < base::min(X1s) ) )
				if( i < j )
				{
					cdfY1[(j-i+1):j] = cdfY0[1:i]
					cdfY1[1:(j-i+1)] = 0.
				}
				else
				{
					cdfY1[1:j] = cdfY0[(i-j+1):i]
				}
			}
			
			if( cdfY1[base::length(cdfY1)] < 1 )
			{
				i = base::max( base::which( x < rvY0$icdf(cdfY1[length(cdfY1)]) ) )
				j = base::min( base::which( cdfY1[1:(length(cdfY1)-1)] == cdfY1[length(cdfY1)] ) )
				if( j > 0 )
				{
					dif = min( length(x) - j , length(x) - i )
					cdfY1[j:(j+dif)] = cdfY0[i:(i+dif)]
					if( j + dif < length(x) ) cdfY1[(j+dif):length(cdfY1)] = 1.
				}
				else
				{
					dev_ok = FALSE
					xdev = 2 * xdev
				}
			}
		}
		
		cdfY1_fct = approxfun( cdfY1 , x )
		Y1 = cdfY1_fct( stats::runif( n = 10000 , min = base::min(cdfY1) , max = base::max(cdfY1) ) )
		
		return(Y1)
	}
	##}}}
	
	)
	##}}}
)

