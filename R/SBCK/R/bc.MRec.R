
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

#' MRec (Matrix Recorrelation) method
#'
#' @description
#' Perform a multivariate bias correction with Gaussian assumption.
#'
#' @details
#' Only pearson correlations are corrected.
#'
#' @references Bárdossy, A. and Pegram, G.: Multiscale spatial recorrelation of
#'             RCM precipitation to produce unbiased climate change scenarios
#'             over large areas and small, Water Resources Research, 48, 9502–,
#'             https://doi.org/10.1029/2011WR011524, 2012.
#'
#' @examples
#' ## Three bivariate random variables (rnorm and rexp are inverted between ref
#' ## and bias)
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
	
	#' @field n_features [integer] Numbers of features
	n_features = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ## {{{
	#' @description
    #' Create a new MRec object.
	#' @param distX [A list of ROOPSD distribution or NULL] Describe the law of
	#'        each margins. A list permit to use different laws for each
	#'        margins. Default is empirical.
	#' @param distY [A list of ROOPSD distribution or NULL] Describe the law of
	#'        each margins. A list permit to use different laws for each
	#'        margins. Default is empirical.
	#'
	#' @return A new `MRec` object.
	initialize = function( distY = NULL , distX = NULL )
	{
		private$distY = distY
		private$distX = distX
	},
	##}}}
	
	## fit ## {{{
	#' @description
    #' Fit the bias correction method
    #' @param Y0 [matrix: n_samples * n_features] Observations in calibration
    #' @param X0 [matrix: n_samples * n_features] Model in calibration
    #' @param X1 [matrix: n_samples * n_features] Model in projection
    #'
    #' @return NULL
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
				private$distX[[i]] = ROOPSD::rv_histogram
		}
		if( is.null(private$distY) )
		{
			private$distY = list()
			for( i in 1:self$n_features)
				private$distY[[i]] = ROOPSD::rv_histogram
		}
		
		## Goto Gaussian world
		private$qmX0 = QM$new( distX0 = private$distX , distY0 = ROOPSD::Normal$new( mean = 0 , sd = 1 ) )
		private$qmX0$fit( X0 = X0 )
		private$qmX1 = QM$new( distX0 = private$distX , distY0 = ROOPSD::Normal$new( mean = 0 , sd = 1 ) )
		private$qmX1$fit( X0 = X1 )
		private$qmY0 = QM$new( distX0 = private$distY , distY0 = ROOPSD::Normal$new( mean = 0 , sd = 1 ) )
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
		private$qmY0 = QM$new( distX0 = ROOPSD::Normal , distY0 = private$distY )
		private$qmY0$fit( Y0 , X0_recor )
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
	##}}}
	
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
#		distX[[i]] = if(!ratio[i]) ROOPSD::Empirical else ROOPSD::EmpiricalRatio
#	
#	
#	## Goto Gaussian world
#	qmX0 = QM$new( distX = distX , distY = ROOPSD::Normal$new( mean = 0 , sd = 1 ) )
#	qmX0$fit( X0 = X0 )
#	qmX1 = QM$new( distX = distX , distY = ROOPSD::Normal$new( mean = 0 , sd = 1 ) )
#	qmX1$fit( X0 = X1 )
#	qmY0 = QM$new( distX = distX , distY = ROOPSD::Normal$new( mean = 0 , sd = 1 ) )
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
#	qmX0Y0 = QM$new( distX = ROOPSD::Normal , distY = distX )
#	qmX0Y0$fit( Y0 , X0_recor )
#	
#	Z0 = qmX0Y0$predict(X0_recor)
#	Z1 = qmX0Y0$predict(X1_recor)
#	
#	return( list( Z0 = Z0 , Z1 = Z1 ) )
#}
#


