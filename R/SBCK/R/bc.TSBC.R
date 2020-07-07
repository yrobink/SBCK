
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

#' TSBC (Time Shifted Bias Correction)
#'
#' Perform a bias correction of auto-correlation
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param lag [integer]
#'        max lag of autocorrelation
#' @param bc_method [SBCK::BC_METHOD]
#'        bias correction method to use after shift of data, default is OTC
#' @param method [character]
#'        If inverse is by row or column, see class Shift
#' @param ref [integer]
#'        reference column/row to inverse shift, see class Shift. Default is 0.5 * (lag+1)
#' @param ...
#'        All others arguments are passed to bc_method
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
#'   \item{\code{new(lag,bc_method,method,ref,...)}}{This method is used to create object of this class with \code{TSBC}}
#'   \item{\code{fit(Y0,X0)}}{Fit the bias correction model X0 with respect to Y0}.
#'   \item{\code{predict(X0)}}{Perform the bias correction of X0.}.
#' }
#' @examples
#' 
#' ## Generate two AR processes
#' X0 = as.vector( stats::arima.sim( n = 2000 , model = list( ar = base::c(  0.6 , 0.2 , -0.1 ) ) , rand.gen = function(n) { return(stats::rnorm( n , mean = 0.2 , sd = 1   )) } )     )
#' Y0 = as.vector( stats::arima.sim( n = 2000 , model = list( ar = base::c( -0.3 , 0.4 , -0.2 ) ) , rand.gen = function(n) { return(stats::rnorm( n , mean = 0   , sd = 0.7 )) } ) + 5 )
#' 
#' ## And correct it with 30 lags
#' tsbc = SBCK::TSBC$new( 30 )
#' tsbc$fit( Y0 , X0 )
#' Z0 = tsbc$predict(X0)
#'
#' @export
TSBC = R6::R6Class( "TSBC" ,
	
	active = list(
	
	method = function(value)
	{
		if(missing(value))
			return(self$shift$method)
		else
			self$shift$method = value
	},
	
	ref = function(value)
	{
		if(missing(value))
			return(self$shift$ref)
		else
			self$shift$ref = value
	}
	
	),
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	shift     = NULL,
	bc_method = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( lag , bc_method = SBCK::OTC , method = "row" , ref = "middle" , ... )
	{
		bc_method_args = list(...)
		self$bc_method = base::do.call( bc_method$new , bc_method_args )
		if( ref == "middle" )
			ref = as.integer(0.5 * (lag+1) )
		self$shift = SBCK::Shift$new( lag , method , ref )
	},
	
	fit = function( Y0 , X0 )
	{
		Xs = self$shift$transform(X0)
		Ys = self$shift$transform(Y0)
		self$bc_method$fit( Ys , Xs )
	},
	
	predict = function(X0)
	{
		return(self$shift$inverse( self$bc_method$predict( self$shift$transform(X0) ) )) 
	}
	
	
	)
)
