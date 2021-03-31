
## Copyright(c) 2021 Yoann Robin
## 
## This file is part of SBCK.
## 
## SBCK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## SBCK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with SBCK.  If not, see <https://www.gnu.org/licenses/>.

#' MBCn (Multivariate Bias Correction)
#'
#' @description
#' Perform a multivariate bias correction.
#'
#' @details
#' BC is performed with an alternance of rotation and univariate BC.
#'
#' @references Cannon, A. J., Sobie, S. R., and Murdock, T. Q.: Bias correction
#'             of simulated precipitation by quantile mapping: how well do
#'             methods preserve relative changes in quantiles and extremes?, J.
#'             Climate, 28, 6938–6959,
#'             https://doi.org/10.1175/JCLI-D-14- 00754.1, 2015.
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
#' XY = SBCK::dataset_gaussian_exp_2d(200)
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
	
	#' @field n_features [integer] Numbers of features
	n_features = NULL,
	#' @field bc [BC class] Univariate BC method
	bc         = NULL,
	#' @field metric [function] distance between two datasets
	metric     = NULL,
	#' @field iter_slope [Stopping class criteria] class used to test when stop
	iter_slope = NULL,
	#' @field bc_params [list] Parameters of bc
	bc_params  = NULL,
	#' @field ortho_mat [array] Array of orthogonal matrix
	ortho_mat  = NULL,
	#' @field tips [array] Array which contains the product of ortho and inverse
	#'        of next
	tips       = NULL,
	#' @field lbc [list] list of BC method used.
	lbc        = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new MBCn object.
	#' @param bc [BC class] Univariate bias correction method
	#' @param metric [function] distance between two datasets
	#' @param stopping_criteria [Stopping class criteria] class use to test when
	#'        to stop the iterations
	#' @param stopping_criteria_params [list] parameters passed to
	#'        stopping_criteria class
	#' @param ... [] Others arguments passed to bc.
	#'
	#' @return A new `MBCn` object.
	initialize = function( bc = QDM , metric = wasserstein , stopping_criteria = SlopeStoppingCriteria , stopping_criteria_params = list( minit = 20 , maxit = 100 , tol = 1e-3 ) , ... ) 
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
	
	## fit ##{{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #'
    #' @return NULL
	fit = function( Y0 , X0 , X1 )
	{
		if( !is.matrix(Y0) ) Y0 = base::matrix( Y0 , ncol = 1 , nrow = length(Y0) )
		if( !is.matrix(X0) ) X0 = base::matrix( X0 , ncol = 1 , nrow = length(X0) )
		if( !is.matrix(X1) ) X1 = base::matrix( X1 , ncol = 1 , nrow = length(X1) )
		self$n_features = base::ncol(Y0)
		
		self$iter_slope$reset()
		maxit = self$iter_slope$maxit
		
		## Generate orthogonal matrix
		self$ortho_mat = ROOPSD::rorthogonal_group( self$n_features , maxit )
		
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
	
	## predict ##{{{
	#' @description
    #' Predict the correction
    #' @param X0 [matrix: n_samples * n_features or NULL] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #'
    #' @return [matrix or list] Return the matrix of correction of X1 if X0 is
    #'                          NULL, else return a list containing Z1 and Z0,
    #'                          the corrections of X1 and X0
	predict = function( X1 , X0 = NULL ) 
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


