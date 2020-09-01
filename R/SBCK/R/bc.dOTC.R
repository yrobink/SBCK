
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

#' dOTC (dynamical Optimal Transport Correction) method
#'
#' @description
#' Perform a multivariate (non stationary) bias correction.
#'
#' @details
#' Three random variables are needed, Y0, X0 and X1. The dynamic between X0 and
#' X1 is estimated, and applied to Y0 to estimate Y1. Finally, OTC is used
#' between X1 and the Y1 estimated.
#'
#' @references Robin, Y., Vrac, M., Naveau, P., Yiou, P.: Multivariate
#'             stochastic bias corrections with optimal transport, Hydrol. Earth
#'             Syst. Sci., 23, 773–786, 2019,
#'             https://doi.org/10.5194/hess-23-773-2019
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#' ## Bin length
#' bin_width = c(0.2,0.2)
#'
#' ## Bias correction
#' ## Step 1 : construction of the class dOTC 
#' dotc = SBCK::dOTC$new( bin_width ) 
#' ## Step 2 : Fit the bias correction model
#' dotc$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction, Z is a list containing
#' ## corrections
#' Z = dotc$predict(X1,X0) 
#' Z$Z0 ## Correction in calibration period
#' Z$Z1 ## Correction in projection period
#'
#' @export
dOTC = R6::R6Class( "dOTC" ,
	inherit = SBCK::OTC,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new dOTC object.
	#' @param bin_width [vector or NULL] A vector of lengths of the cells
	#'        discretizing R^{numbers of variables}. If NULL, it is estimating
	#'        during the fit
	#' @param bin_origin [vector or NULL] Coordinate of lower corner of one
	#'        cell. If NULL, c(0,...,0) is used
	#' @param cov_factor [string or matrix] Covariance factor to correct the
	#'        dynamic transferred between X0 and Y0. For string, available
	#'        values are "std" and "cholesky"
	#' @param ot [OTSolver] Optimal Transport solver, default is the network
	#'        simplex
	#'
	#' @return A new `dOTC` object.
	initialize = function( bin_width = NULL , bin_origin = NULL , cov_factor = "std" , ot = SBCK::OTNetworkSimplex$new() )
	{
		super$initialize( bin_width , bin_origin , ot )
		private$cov_factor = cov_factor
	},
	##}}}
	
	## Fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return NULL
	fit = function( Y0 , X0 , X1 )
	{
		## Dimension and data formating
		if( !is.matrix(Y0) ) Y0 = matrix( Y0 , nrow = length(Y0) , ncol = 1 )
		if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		if( !is.matrix(X1) ) X1 = matrix( X1 , nrow = length(X1) , ncol = 1 )
		
		## Bin width
		if( is.null(self$bin_width) )
		{
			self$bin_width = SBCK::bin_width_estimator( list(Y0,X0,X1) )
		}
		if( is.null(self$bin_origin) )
		{
			self$bin_origin = base::rep( 0. , length(self$bin_width) )
		}
		
		self$n_features = base::ncol(Y0)
		
		## Correction factor of motion
		cf = NULL
		if( self$n_features == 1 && !is.numeric(private$cov_factor) )
		{
			cf = sd(Y0) * sd(X0)^(-1)
		}
		else
		{
			if( !is.matrix(private$cov_factor) )
			{
				if( private$cov_factor == "cholesky" )
				{
					cf = base::chol( stats::cov(Y0) ) %*% base::solve( base::chol( stats::cov(X0) ) )
				}
				else
				{
					cf = base::diag( base::apply( Y0 , 2 , sd ) * base::apply( X0 , 2 , stats::sd )^(-1) )
				}
			}
		}
		## Motion
		private$otcB0R0 = SBCK::OTC$new( self$bin_width , self$bin_origin ) ## Bias
		otcB0B1 = SBCK::OTC$new( self$bin_width , self$bin_origin ) ## Evolution
		private$otcB0R0$fit( Y0 , X0 )
		otcB0B1$fit( X1 , X0 )
		
		motion = otcB0B1$predict(X0) - X0
		if( self$n_features == 1 )
		{
			motion = cf * motion
		}
		else
		{
			motion = base::t(base::apply( motion , 1 , function(x){ invisible( cf %*% x ) } ))
		}
		
		## Estimation of Y1
		Y1 = private$otcB0R0$predict( X0 ) + motion
		
		super$fit( Y1 , X1 )
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #'
	#' Note: Only the center of the bins associated to the corrected points are
	#' returned, but all corrections of the form:
	#' >> bw = dotc$bin_width / 2
	#' >> n  = base::prod(base::dim(X1))
	#' >> Z1 = dotc$predict(X1)
	#' >> Z1 = Z1 + t(matrix(stats::runif( n = n min = - bw , max = bw ) , ncol = dim(X1)[1] ))
	#' are equivalent for OTC.
	#'
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 , X0 = NULL )
	{
		Z1 = super$predict(X1)
		
		if( !is.null(X0) )
		{
			Z0 = private$otcB0R0$predict(X0)
			return( list( Z1 = Z1 , Z0 = Z0 ) )
		}
		return(Z1)
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
	
	otcB0R0 = NULL ,
	cov_factor = NULL
	)
)
