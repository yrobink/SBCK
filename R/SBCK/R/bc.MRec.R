
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

#' MRec (Matrix Recorrelation) method
#'
#' Perform a multivariate bias correction with Gaussian assumption. Only pearson correlations are corrected.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param distX [A list of ROOPSD_ distribution or NULL]
#'        Describe the law of each margins. A list permit to use different laws for each margins. Default is empirical.
#' @param distY [A list of ROOPSD_ distribution or NULL]
#'        Describe the law of each margins. A list permit to use different laws for each margins. Default is empirical.
#' @param Y0  [matrix]
#'        A matrix containing references during calibration period (time in column, variables in row)
#' @param X0 [matrix]
#'        A matrix containing biased data during calibration period (time in column, variables in row)
#' @param X1 [matrix]
#'        A matrix containing biased data during projection period (time in column, variables in row)
#'
#' @return Object of \code{\link{R6Class}} with methods for bias correction
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(distY,distX)}}{This method is used to create object of this class with \code{MRec}}
#'   \item{\code{fit(Y0,X0,X1)}}{Fit the bias correction model from Y0, X0 and X1}.
#' }
#' @references Bárdossy, A. and Pegram, G.: Multiscale spatial recorrelation of RCM precipitation to produce unbiased climate change scenarios over large areas and small, Water Resources Research, 48, 9502–, https://doi.org/10.1029/2011WR011524, 2012.
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#' ## Bias correction
#' ## Step 1 : construction of the class MRec 
#' mrec = SBCK::MRec$new() 
#' ## Step 2 : Fit the bias correction model
#' mrec$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction, Z is a list containing corrections.
#' Z = mrec$predict(X1,X0) ## X0 is optional, in this case Z0 is NULL
#' Z$Z0 ## Correction in calibration period
#' Z$Z1 ## Correction in projection period
#'
#' @export
MRec = R6::R6Class( "MRec" ,
	
	
	public = list(
	
	
	###############
	## Arguments ##
	###############
	
	n_features = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( distY = NULL , distX = NULL )
	{
		private$distY = distY
		private$distX = distX
	},
	
	fit = function( Y0 , X0 , X1 )
	{
		## Data in matrix
		if( !is.matrix(Y0) ) Y0 = base::matrix( Y0 , ncol = 1 , nrow = length(Y0) )
		if( !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		if( !is.matrix(X1) ) X1 = base::matrix( X1 , ncol = 1 , nrow = length(X1) )
		self$n_features = base::ncol(Y0)
		
		## Kind of variable
		if( is.null(private$distX) )
		{
			private$distX = list()
			for( i in 1:self$n_features)
				private$distX[[i]] = ROOPSD_rv_histogram
		}
		if( is.null(private$distY) )
		{
			private$distY = list()
			for( i in 1:self$n_features)
				private$distY[[i]] = ROOPSD_rv_histogram
		}
		
		## Goto Gaussian world
		private$qmX0 = QM$new( distX0 = private$distX , distY0 = ROOPSD_Normal$new( mean = 0 , sd = 1 ) )
		private$qmX0$fit( X0 = X0 )
		private$qmX1 = QM$new( distX0 = private$distX , distY0 = ROOPSD_Normal$new( mean = 0 , sd = 1 ) )
		private$qmX1$fit( X0 = X1 )
		private$qmY0 = QM$new( distX0 = private$distY , distY0 = ROOPSD_Normal$new( mean = 0 , sd = 1 ) )
		private$qmY0$fit( X0 = Y0 )
		Y0g = private$qmY0$predict(Y0)
		X0g = private$qmX0$predict(X0)
		X1g = private$qmX1$predict(X1)
		
		## Correlation
		CY0g = stats::cor( Y0g , method = "pearson" )
		CX0g = stats::cor( X0g , method = "pearson" )
		
		## Squareroot
		svdY0g  = base::svd(CY0g)
		private$S_CY0g  = svdY0g$u %*% base::diag(base::sqrt(svdY0g$d)) %*% base::t(svdY0g$u)
		svdX0g  = base::svd(CX0g)
		private$Si_CX0g = svdX0g$u %*% base::diag(1./base::sqrt(svdX0g$d)) %*% base::t(svdX0g$u)
		private$re_un_mat = private$S_CY0g %*% private$Si_CX0g
		
		## Decor-recor-relation
		X0_recor = base::t(private$re_un_mat %*% base::t(X0g))
		X1_recor = base::t(private$re_un_mat %*% base::t(X1g))
		
		## Final QM
		private$qmY0 = QM$new( distX0 = ROOPSD_Normal , distY0 = private$distY )
		private$qmY0$fit( Y0 , X0_recor )
	},
	
	predict = function( X1 , X0 = NULL )
	{
		X1g = private$qmX1$predict(X1)
		X1_recor = base::t(private$re_un_mat %*% base::t(X1g))
		Z1 = private$qmY0$predict(X1_recor)
		
		Z0 = NULL
		if( !is.null(X0) )
		{
			X0g = private$qmX0$predict(X0)
			X0_recor = base::t(private$re_un_mat %*% base::t(X0g))
			Z0 = private$qmY0$predict(X0_recor)
			return( list( Z1 = Z1 , Z0 = Z0 ) )
		}
		return(Z1)
	}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	distX = NULL,
	distY = NULL,
	S_CY0g  = NULL,
	Si_CX0g = NULL,
	re_un_mat = NULL,
	qmX0 = NULL,
	qmX1 = NULL,
	qmY0 = NULL
	
	)
)



#MRec = function( X0 , X1 , Y0 , ratio = NULL )
#{
#	## Data in matrix
#	if( !is.matrix(Y0) ) Y0 = base::matrix( Y0 , ncol = 1 , nrow = length(Y0) )
#	if( !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
#	if( !is.matrix(X1) ) X1 = base::matrix( X1 , ncol = 1 , nrow = length(X1) )
#	n_features = base::ncol(Y0)
#	
#	## Kind of variable
#	if( is.null(ratio) )
#		ratio = base::rep( FALSE , n_features )
#	
#	distX = list()
#	for( i in 1:n_features)
#		distX[[i]] = if(!ratio[i]) ROOPSD_Empirical else ROOPSD_EmpiricalRatio
#	
#	
#	## Goto Gaussian world
#	qmX0 = QM$new( distX = distX , distY = ROOPSD_Normal$new( mean = 0 , sd = 1 ) )
#	qmX0$fit( X0 = X0 )
#	qmX1 = QM$new( distX = distX , distY = ROOPSD_Normal$new( mean = 0 , sd = 1 ) )
#	qmX1$fit( X0 = X1 )
#	qmY0 = QM$new( distX = distX , distY = ROOPSD_Normal$new( mean = 0 , sd = 1 ) )
#	qmY0$fit( X0 = Y0 )
#	Y0g = qmY0$predict(Y0)
#	X0g = qmX0$predict(X0)
#	X1g = qmX1$predict(X1)
#	
#	## Correlation
#	CY0g = stats::cor( Y0g , method = "pearson" )
#	CX0g = stats::cor( X0g , method = "pearson" )
#	
#	## Squareroot
#	svdY0g  = base::svd(CY0g)
#	S_CY0g  = svdY0g$u %*% base::diag(base::sqrt(svdY0g$d)) %*% base::t(svdY0g$u)
#	svdX0g  = base::svd(CX0g)
#	Si_CX0g = svdX0g$u %*% base::diag(1./base::sqrt(svdX0g$d)) %*% base::t(svdX0g$u)
#	
#	## Decor-recor-relation
#	X0_recor = base::t(S_CY0g %*% Si_CX0g %*% base::t(X0g))
#	X1_recor = base::t(S_CY0g %*% Si_CX0g %*% base::t(X1g))
#	
#	## Final QM
#	qmX0Y0 = QM$new( distX = ROOPSD_Normal , distY = distX )
#	qmX0Y0$fit( Y0 , X0_recor )
#	
#	Z0 = qmX0Y0$predict(X0_recor)
#	Z1 = qmX0Y0$predict(X1_recor)
#	
#	return( list( Z0 = Z0 , Z1 = Z1 ) )
#}
#


