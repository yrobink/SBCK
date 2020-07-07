
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
#' Perform an univariate bias correction of X0 with respect to Y0 (correction is applied margins by margins).
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param distX0 [A ROOPSD_ distribution or a list of them]
#'        Describe the law of each margins. A list permit to use different laws for each margins. Default is empirical.
#' @param distY0 [A ROOPSD_ distribution or a list of them]
#'        Describe the law of each margins. A list permit to use different laws for each margins. Default is empirical.
#' @param ...
#'        Others optional named arguments:
#' @param n_features  [integer]
#'        Normaly infered during fit, but if distX0 and distY are simultaneously frozen, must be set during initialization.
#' @param kwargsX0  [list]
#'        Arguments passed to distX0
#' @param kwargsY0  [list]
#'        Arguments passed to distY0
#' @param Y0  [matrix]
#'        A matrix containing references during calibration period (time in column, variables in row)
#' @param X0 [matrix]
#'        A matrix containing biased data during calibration period (time in column, variables in row)
#'
#' @return Object of \code{\link{R6Class}} with methods for bias correction
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(distX0,distY0,...)}}{This method is used to create object of this class with \code{QM}}
#'   \item{\code{fit(Y0,X0)}}{Fit the bias correction model from Y0 and X0}.
#'   \item{\code{predict(X0)}}{Perform the bias correction of X0 with respect to Y0.}.
#' }
#' @references Panofsky, H. A. and Brier, G. W.: Some applications of statistics to meteorology, Mineral Industries Extension Services, College of Mineral Industries, Pennsylvania State University, 103 pp., 1958.
#' @references Wood, A. W., Leung, L. R., Sridhar, V., and Lettenmaier, D. P.: Hydrologic Implications of Dynamical and Statistical Approaches to Downscaling Climate Model Outputs, Clim. Change, 62, 189–216, https://doi.org/10.1023/B:CLIM.0000013685.99609.9e, 2004.
#' @references Déqué, M.: Frequency of precipitation and temperature extremes over France in an anthropogenic scenario: Model results and statistical correction according to observed values, Global Planet. Change, 57, 16–26, https://doi.org/10.1016/j.gloplacha.2006.11.030, 2007.
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#'
#' ## Bias correction
#' ## Step 1 : construction of the class QM 
#' qm = SBCK::QM$new() 
#' ## Step 2 : Fit the bias correction model
#' qm$fit( Y0 , X0 )
#' ## Step 3 : perform the bias correction, Z0 is the correction of
#' ## X0 with respect to the estimation of Y0
#' Z0 = qm$predict(X0)
#'
#' # ## But in fact the laws are known, we can fit parameters:
#' distY0 = list( ROOPSD_Exponential , ROOPSD_Normal )
#' distX0 = list( ROOPSD_Normal , ROOPSD_Exponential )
#' qm_fix = SBCK::QM$new( distY0 = distY0 , distX0 = distX0 )
#' qm_fix$fit( Y0 , X0 )
#' Z0 = qm_fix$predict(X0) 
#' @export
QM = R6::R6Class( "QM",
	## Public elements
	##============={{{
	
	public = list(
	
	## Arguments
	##==========
	distX0 = NULL,
	distY0 = NULL,
	n_features = NULL,
	tol = 1e-3,
	
	## Constructor
	##============
	initialize = function( distX0 = ROOPSD_rv_histogram , distY0 = ROOPSD_rv_histogram , ... )##{{{
	{
		kwargs = list(...)
		self$distX0 = DistHelper$new( dist = distX0 , kwargs = kwargs[["kwargsX0"]] )
		self$distY0 = DistHelper$new( dist = distY0 , kwargs = kwargs[["kwargsY0"]] )
		self$n_features = kwargs[["n_features"]]
		self$tol = if( is.null(kwargs[["tol"]]) ) 1e-3 else kwargs[["tol"]]
	},
	##}}}
	
	
	## Methods
	##========
	
	fit = function( Y0 = NULL , X0 = NULL ) ##{{{
	{
		## Data in matrix
		if( !is.null(Y0) && !is.matrix(Y0) ) Y0 = base::matrix( Y0 , ncol = 1 , nrow = length(Y0) )
		if( !is.null(X0) && !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		
		## Find n_features
		if( is.null(self$n_features ) )
		{
			if( !is.null(Y0) )
			{
				self$n_features = base::ncol(Y0)
			}
			else if( !is.null(X0) )
			{
				self$n_features = base::ncol(X0)
			}
			else
			{
				base::stop( "QM fit: if X0 = Y0 = NULL, n_features must be set during intialization" )
			}
		}
		
		## distX and distY into list
		self$distY0$set_features(self$n_features)
		self$distX0$set_features(self$n_features)
		
		## Now fit itself
		for( i in 1:self$n_features )
		{
			self$distY0$fit( Y0[,i] , i )
			self$distX0$fit( X0[,i] , i )
		}
		
	},
	##}}}
	
	predict = function(X0)##{{{
	{
		if( !is.null(X0) && !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		Z0 = base::matrix( NA , nrow = base::nrow(X0) , ncol = base::ncol(X0) )
		for( i in 1:self$n_features )
		{
			cdf = self$distX0$law[[i]]$cdf( X0[,i] )
			cdf[!(cdf < 1)] = 1-self$tol
			cdf[!(cdf > 0)] = self$tol
			Z0[,i] = self$distY0$law[[i]]$icdf( cdf )
		}
		return(Z0)
	}
	##}}}
	
	),
	##}}}
	
	## Private elements
	##==============={{{
	private = list(
	)
	##}}}

)

