
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

#' MBCn (Multivariate Bias Correction)
#'
#' Perform a multivariate bias correction.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param bc [bias correction methd]
#'        Non stationary BC method of SBCK, as QDM. Default is QDM
#' @param metric [function]
#'        Distance between two distributions. Default is wasserstein.
#' @param stopping_criteria [R6]
#'        Class which implement a criteria to stop iterations. See SlopeStoppingCriteria
#' @param stopping_criteria_params [list]
#'        Params
#' @param ... 
#'        Named arguments passed to bc method
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
#'   \item{\code{new(bin_width,bin_origin,cov_factor)}}{This method is used to create object of this class with \code{MBCn}}
#'   \item{\code{fit(Y0,X0,X1)}}{Fit the bias correction model from Y0, X0 and X1}.
#'   \item{\code{predict(X1,X0)}}{Perform the bias correction.}
#' }
#' @references Cannon, A. J., Sobie, S. R., and Murdock, T. Q.: Bias correction of simulated precipitation by quantile mapping: how well do methods preserve relative changes in quantiles and extremes?, J. Climate, 28, 6938–6959, https://doi.org/10.1175/JCLI-D-14- 00754.1, 2015.
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(2000)
#' X0 = XY$X0 ## Biased in calibration period
#' Y0 = XY$Y0 ## Reference in calibration period
#' X1 = XY$X1 ## Biased in projection period
#'
#' ## Bias correction
#' ## Step 1 : construction of the class MBCn
#' mbcn = SBCK::MBCn$new() 
#' ## Step 2 : Fit the bias correction model
#' mbcn$fit( Y0 , X0 , X1 )
#' ## Step 3 : perform the bias correction, Z is a list containing
#' ## corrections
#' Z = mbcn$predict(X1,X0) 
#' Z$Z0 ## Correction in calibration period
#' Z$Z1 ## Correction in projection period
#'
#' @export
MBCn = R6::R6Class( "MBCn" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	n_features = NULL,
	bc = NULL,
	metric = NULL,
	iter_slope = NULL,
	bc_params = NULL,
	ortho_mat = NULL,
	tips = NULL,
	lbc = NULL,
	
	#################
	## Constructor ##
	#################
	
	initialize = function( bc = QDM , metric = wasserstein , stopping_criteria = SlopeStoppingCriteria , stopping_criteria_params = list( minit = 20 , maxit = 100 , tol = 1e-3 ) , ... ) ##{{{
	{
		self$n_features = NULL
		self$bc = bc
		self$metric = metric
		self$iter_slope = base::do.call( stopping_criteria$new , stopping_criteria_params )
		self$bc_params = list(...)
		self$ortho_mat = NULL
		self$tips = NULL
		self$lbc = list()
	},
	##}}}
	
	fit = function( Y0 , X0 , X1 )##{{{
	{
		if( !is.matrix(Y0) ) Y0 = base::matrix( Y0 , ncol = 1 , nrow = length(Y0) )
		if( !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		if( !is.matrix(X1) ) X1 = base::matrix( X1 , ncol = 1 , nrow = length(X1) )
		self$n_features = base::ncol(Y0)
		
		self$iter_slope$reset()
		maxit = self$iter_slope$maxit
		
		## Generate orthogonal matrix
		self$ortho_mat = rorthogonal_group( self$n_features , maxit )
		
		## Tips for performance: inverse + ortho of next in one pass
		self$tips = array( NA , dim = base::dim(self$ortho_mat) )
		for( i in 1:(maxit-1) )
			self$tips[,,i] = self$ortho_mat[,,i+1] %*% base::solve(self$ortho_mat[,,i])
		self$tips[,,maxit] = base::solve(self$ortho_mat[,,maxit])
		
		## Init loop
		Z0_o = base::t(self$ortho_mat[,,1] %*% base::t(X0))
		Z1_o = base::t(self$ortho_mat[,,1] %*% base::t(X1))
		
		## Main loop
		while(!self$iter_slope$stop)
		{
			nit = self$iter_slope$nit
			Y0_o = base::t(self$ortho_mat[,,nit] %*% base::t(Y0))
			
			bc = base::do.call( self$bc$new , self$bc_params )
			bc$fit( Y0_o , Z0_o , Z1_o )
			Z = bc$predict(Z1_o,Z0_o)
			Z1_o = Z$Z1
			Z0_o = Z$Z0
			self$lbc[[nit]] = bc
			
			self$iter_slope$append(self$metric(Z0_o,Y0_o))
			
			Z0_o = base::t(self$tips[,,nit] %*% base::t(Z0_o))
			Z1_o = base::t(self$tips[,,nit] %*% base::t(Z1_o))
		}
		
		nit = self$iter_slope$nit
		self$ortho_mat = self$ortho_mat[,,1:nit]
		self$tips = self$tips[,,1:nit]
		self$tips[,,nit] = base::solve(self$ortho_mat[,,nit]) 
		
		Z0 = base::t(self$tips[,,nit] %*% base::t(Z0_o))
		Z1 = base::t(self$tips[,,nit] %*% base::t(Z1_o))
		
		bc = base::do.call( self$bc$new , self$bc_params )
		bc$fit( Y0 , Z0 , Z1 )
		self$lbc[[nit]] = bc
	},
	##}}}
	
	predict = function( X1 , X0 = NULL ) ##{{{
	{
		if( is.null(X0) )
			return(private$predict_X1(X1))
		else
			return(private$predict_X1_X0(X1,X0))
	
	}
	##}}}
	
	),
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	
	#############
	## Methods ##
	#############
	
	predict_X1 = function(X1) ##{{{
	{
		if( !is.matrix(X1) ) X1 = base::matrix( X1 , ncol = 1 , nrow = length(X1) )
		
		nit = self$iter_slope$nit
		
		Z1_o = base::t(self$ortho_mat[,,1] %*% base::t(X1))
		
		for( i in 1:(nit-1) )
		{
			Z1_o = self$lbc[[i]]$predict(Z1_o)
			Z1_o = base::t(self$tips[,,i] %*% base::t(Z1_o))
		}
		
		Z1_o = base::t(self$tips[,,nit] %*% base::t(Z1_o))
		Z1 = self$lbc[[nit]]$predict(Z1_o)
		return(Z1)
		
	},
	##}}}
	
	predict_X1_X0 = function(X1,X0) ##{{{
	{
		if( !is.matrix(X1) ) X1 = base::matrix( X1 , ncol = 1 , nrow = length(X1) )
		if( !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		
		nit = self$iter_slope$nit
		
		Z0_o = base::t(self$ortho_mat[,,1] %*% base::t(X0))
		Z1_o = base::t(self$ortho_mat[,,1] %*% base::t(X1))
		
		for( i in 1:(nit-1) )
		{
			Z = self$lbc[[i]]$predict(Z1_o,Z0_o)
			Z0_o = base::t(self$tips[,,i] %*% base::t(Z$Z0))
			Z1_o = base::t(self$tips[,,i] %*% base::t(Z$Z1))
		}
		
		Z0_o = base::t(self$tips[,,nit] %*% base::t(Z0_o))
		Z1_o = base::t(self$tips[,,nit] %*% base::t(Z1_o))
		Z = self$lbc[[nit]]$predict(Z1_o,Z0_o)
		return(Z)
	}
	##}}}
	
	)
)


