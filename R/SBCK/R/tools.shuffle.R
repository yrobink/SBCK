
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

## SchaakeShuffle ##{{{


#' ShaakeShuffle class
#'
#' @description
#' Perform the Schaake Shuffle 
#'
#' @details
#' as fit/predict mode
#'
#'
#' @examples
#' X0 = matrix( stats::runif(20) , ncol = 2 )
#' Y0 = matrix( stats::runif(20) , ncol = 2 )
#' ss = SchaakeShuffle$new()
#' ss$fit(Y0)
#' Z0 = ss$predict(X0)
#'
#' @export
SchaakeShuffle = R6::R6Class( "SchaakeShuffle" ,
	
	public = list( ##{{{
	
	## initialize ##{{{
	#' @description
    #' Create a new ShaakeShuffle object.
	#' @param Y0 [vector] The reference vector
	#'
	#' @return A new `ShaaleShuffle` object.
	initialize = function( Y0 = NULL )
	{
		private$Y0 = NULL
		if( !is.null(Y0) )
			self$fit(Y0)
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the model
	#' @param Y0 [vector] The reference vector
	#'
	#' @return NULL
	fit = function( Y0 )
	{
		private$Y0 = Y0
		if( !is.matrix(Y0) ) private$Y0 = matrix( Y0 , nrow = length(Y0) , ncol = 1 )
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Fit the model
	#' @param X0 [vector] The vector to apply shuffle
	#'
	#' @return Z0 [vector] data shuffled
	predict = function( X0 )
	{
		if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		
		YY = NULL
		XX = NULL
		nrowY0 = base::nrow(private$Y0)
		nrowX0 = base::nrow(X0)
		n_features = base::ncol(private$Y0)
		if( nrowY0 < nrowX0 )
		{
			YY = matrix( NA , nrow = nrowX0 , ncol = n_features )
			YY[1:nrowY0,] = private$Y0
			YY[nrowY0:nrowX0,] = private$Y0[base::sample( nrowY0 , nrowX0 - nrowY0 , replace = TRUE ),]
			XX = X0
		}
		else if( nrowX0 < nrowY0 )
		{
			XX = matrix( NA , nrow = nrowY0 , ncol = n_features )
			XX[1:nrowX0,] = X0
			XX[nrowX0:nrowY0,] = X0[base::sample( nrowX0 , nrowY0 - nrowX0 , replace = TRUE ),]
			YY = private$Y0
		}
		else
		{
			XX = X0
			YY = private$Y0
		}
		
		ZZ = matrix( NA , nrow = base::nrow(XX) , ncol = base::ncol(XX) )
		for( i in 1:n_features )
		{
			ZZ[,i] = private$predict_univariate( YY[,i] , XX[,i] )
		}
		Z = ZZ[1:nrowX0,]
		
		return(Z)
	}
	##}}}
	
	),
	##}}}
	
	private = list(##{{{
	
	###############
	## Arguments ##
	###############
	
	Y0 = NULL,
	
	
	#############
	## Methods ##
	#############
	
	predict_univariate = function( Y , X )##{{{
	{
		rank_X  = base::rank(X)
		rank_Y  = base::rank(Y)
		arank_X = base::order(rank_X)
		Z = X[arank_X][rank_Y]
		return(Z)
	}
	##}}}
	
	)
	##}}}
)
##}}}

## schaake_shuffle ##{{{

#' schaake_shuffle function
#'
#' Apply the Schaake shuffle to transform the rank of X0 such that its
#' correspond to the rank of Y0
#' 
#' @usage schaake_shuffle(Y0,X0)
#' @param X0 [vector] The vector to transform the rank
#' @param Y0 [vector] The reference vector
#'
#' @return Z0 [vector] X shuffled.
#'
#' @examples
#' X0 = stats::runif(10)
#' Y0 = stats::runif(10)
#' Z0 = SBCK::schaake_shuffle( Y0 , X0 )
#'
#' @export
schaake_shuffle = function( Y0 , X0 )
{
	ss = SchaakeShuffle$new(Y0)
	return(ss$predict(X0))
}
##}}}

## SchaakeShuffleRef ##{{{

#' ShaakeShuffleRef class
#'
#' @description
#' Match the rank structure of X with them of Y by reordering X.
#'
#' @details
#' Fix one features to keep the structure of X.
#'
#' @param ref [integer] The reference
#' @param X0 [vector] The vector to transform the rank
#' @param Y0 [vector] The reference vector
#'
#' @examples
#' X0 = matrix( stats::runif(20) , ncol = 2 )
#' Y0 = matrix( stats::runif(20) , ncol = 2 )
#' ss = SchaakeShuffleRef$new( ref = 1 )
#' ss$fit(Y0)
#' Z0 = ss$predict(X0)
#'
#' @export
SchaakeShuffleRef = R6::R6Class( "SchaakeShuffleRef" ,
	
	inherit = SchaakeShuffle,
	
	public = list( ##{{{
	
	#' @field ref [integer] Reference
	ref = NULL,
	
	## initialize ##{{{
	#' @description
    #' Create a new ShaakeShuffleRef object.
	#' @param ref [integer] Reference
	#' @param Y0 [vector] The reference vector
	#'
	#' @return A new `ShaaleShuffleRef` object.
	initialize = function( ref , Y0 = NULL )
	{
		self$ref = ref
		super$initialize(Y0)
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit the model
	#' @param Y0 [vector] The reference vector
	#'
	#' @return NULL
	fit = function(Y0)
	{
		super$fit(Y0)
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Fit the model
	#' @param X0 [vector] The vector to apply shuffle
	#'
	#' @return Z0 [vector] data shuffled
	predict = function(X0)
	{
		if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		Z0 = super$predict(X0)
		
		rank_X0  = base::rank(X0[,self$ref])
		rank_Z0  = base::rank(Z0[,self$ref])
		arank_Z0 = base::order(rank_Z0)
		Z0 = Z0[arank_Z0,][rank_X0,]
		return(Z0)
	}
	##}}}
	
	)
	##}}}

)
##}}}

## SchaakeShuffleMultiRef ##{{{

#' ShaakeShuffleMultiRef class
#'
#' @description
#' Match the rank structure of X with them of Y by reordering X.
#'
#' @details
#' Can keep multiple features to keep the structure of X.
#'
#' @param lag_search An integer corresponding to the number of time lags to
#' account for when searching in \code{refdata} the best analogue of the
#' conditioning dimension rank association observed in \code{bc1d}. The default
#' value is no lag, i.e., \code{lag_search}=0.
#' @param lag_keep An integer corresponding to the number of time lags to keep
#' in the correction for each best analogue of rank associations found.
#' \code{lag_keep} has to be lower or equal to \code{lag_search}.  The default
#' value is no lag, i.e., \code{lag_search}=0.
#'
#' @examples
#' X0 = matrix( stats::runif(50) , ncol = 2 )
#' Y0 = matrix( stats::runif(50) , ncol = 2 )
#' ssmr = SchaakeShuffleMultiRef$new( lag_search = 3 , lag_keep = 1 , cond_cols = 1 )
#' ssmr$fit(Y0)
#' Z0 = ssmr$predict(X0)
#'
#' @export
SchaakeShuffleMultiRef = R6::R6Class( "SchaakeShuffleMultiRef" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field cond_cols [vector of integer] The conditioning columns
	cond_cols  = NULL,
	#' @field lag_search [integer] Number of lag to take into account
	lag_search = NULL,
	#' @field lag_keep [integer] Number of lag to keep
	lag_keep   = NULL,
	#' @field Y0 [matrix] Reference data
	Y0         = NULL,
	
	#################
	## Constructor ##
	#################
	
	## Init  ##{{{
	#' @description
    #' Create a new ShaakeShuffleMultiRef object.
	#' @param cond_cols [vector of integer] The conditioning columns
	#' @param lag_search [integer] Number of lag to take into account
	#' @param lag_keep [integer] Number of lag to keep
	#'
	#' @return A new `ShaaleShuffleMultiRef` object.
	initialize = function( lag_search , lag_keep , cond_cols = base::c(1) )
	{
		self$cond_cols  = cond_cols
		self$lag_search = lag_search
		self$lag_keep   = lag_keep
	},
	##}}}
	
	## Fit  ##{{{
	#' @description
    #' Fit the model
	#' @param Y0 [vector] The reference vector
	#'
	#' @return NULL
	fit = function(Y0)
	{
		self$Y0 = Y0
		if( !is.matrix(Y0) ) self$Y0 = matrix( Y0 , nrow = length(Y0) , ncol = 1 )
	},
	##}}}
	
	## predict  ##{{{
	#' @description
    #' Fit the model
	#' @param X0 [vector] The vector to apply shuffle
	#'
	#' @return Z0 [vector] data shuffled
	predict = function(X0)
	{
		## Reorganize data
		if( !is.matrix(X0) ) X0 = matrix( X0 , nrow = length(X0) , ncol = 1 )
		
		## If nrow(X0) != nrow(Y0), extend data with a resampling
		YY = NULL
		XX = NULL
		nrowY0 = base::nrow(self$Y0)
		nrowX0 = base::nrow(X0)
		n_features = base::ncol(self$Y0)
		if( nrowY0 < nrowX0 )
		{
			YY = matrix( NA , nrow = nrowX0 , ncol = n_features )
			YY[1:nrowY0,] = private$Y0
			YY[nrowY0:nrowX0,] = self$Y0[base::sample( nrowY0 , nrowX0 - nrowY0 , replace = TRUE ),]
			XX = X0
		}
		else if( nrowX0 < nrowY0 )
		{
			XX = matrix( NA , nrow = nrowY0 , ncol = n_features )
			XX[1:nrowX0,] = X0
			XX[nrowX0:nrowY0,] = X0[base::sample( nrowX0 , nrowY0 - nrowX0 , replace = TRUE ),]
			YY = self$Y0
		}
		else
		{
			XX = X0
			YY = self$Y0
		}
		n_samples  = nrow(XX)
		
		ZZ   = array(NaN, dim = dim(XX)) # (N days x P var)
		ZZ_s = apply(XX, 2, sort)
		ZZ_r = apply(XX, 2, rank, ties.method = "min")
		YY_r = apply(YY, 2, rank, ties.method = "min")
		
		ranks_conddim_REF = YY_r[, self$cond_cols, drop = FALSE]
		ranks_conddim_BC  = ZZ_r[, self$cond_cols, drop = FALSE]
		
		nsearch <- ((n_samples - self$lag_search - 1) %/% (self$lag_keep + 1)) + (((n_samples - self$lag_search - 1) %% (self$lag_keep + 1)) > 0) + 1
		
		time_bestanalogue = numeric(length = nsearch)
		dist_bestanalogue = numeric(length = nsearch)
		visited_time      = numeric(length = n_samples)
		
		## block search in matrix
		iBS = toeplitz(1:n_samples)[seq(self$lag_search+1,n_samples,self$lag_keep+1),(self$lag_search+1):1]
		if( iBS[nrow(iBS),ncol(iBS)] < n_samples )
			iBS = rbind( iBS , (n_samples-self$lag_search):n_samples )
		
		## pairwise distances
		A = toeplitz(as.vector(ranks_conddim_REF))[(ncol(iBS)):nrow(ranks_conddim_REF),(ncol(iBS)):1]
		B = matrix( ranks_conddim_BC[as.vector(iBS),] , ncol = ncol(iBS) )
		D = SBCK::pairwise_distances( B , A )
		
		## Correction
		dist_be  = base::apply( D , 1 , base::min )
		time_be  = base::apply( D , 1 , base::which.min ) + ncol(iBS) - 1
		vis_time = numeric(n_samples)
		
		## The first line
		ibk = iBS[1,]
		ibr = (time_be[1] - length(ibk)+1):time_be[1]
		for( i in 1:n_features )
			ZZ[ibk,i] = ZZ_s[YY_r[ibr,i],i]
		
		## iblockref and iblockkeep in matrix, except the first line
		iBK = iBS[2:nrow(iBS),(ncol(iBS)-self$lag_keep):ncol(iBS)]
		iBR = base::t(apply( matrix(time_be[2:length(time_be)],nrow = length(time_be) -1 ) , 1 , function(x) { (x-self$lag_keep):x } ))
		## RR = Rank Ref
		RR = YY_r[as.vector(base::t(iBR)),]
		## Use rbind to add false index to have the same dimension that ZZ_s, add k*n_samples at each col to have the good uni-dim index
		RR = base::t(base::t(rbind( RR , matrix( 1 , nrow = self$lag_search + 1 , ncol = n_features ) )) + seq(0,n_samples*(n_features-1),n_samples))
		## Now RR can be used in uni-dim, and remove last "false" index.
		ZZ[as.vector(base::t(iBK)),] = matrix( ZZ_s[as.vector(RR)] , ncol = n_features )[-seq(n_samples-self$lag_search,n_samples),]
		
		## Compute visited time
		table_visited = table(as.vector(iBR))
		idx           = as.integer(names(table_visited))
		vis_time[idx] = as.vector(table_visited)
		vis_time[ibr] = vis_time[ibr] + 1
		
		## Shuffled data
		Z = ZZ[1:nrowX0,]
		
		return(Z)
	}
	##}}}
	
	)
	
)
##}}}

## MVQuantilesShuffle ##{{{

#' MVQuantilesShuffle
#'
#' @description
#' Multivariate Schaake shuffle using the quantiles.
#'
#' @details
#' Used to reproduce the dependence structure of a dataset to another dataset
#'
#' @references Vrac, M. et S. Thao (2020). “R2 D2 v2.0 : accounting for temporal
#'             dependences in multivariate bias correction via analogue rank
#'             resampling”. In : Geosci. Model Dev. 13.11, p. 5367-5387.
#'             doi :10.5194/gmd-13-5367-2020.
#'
#' @importFrom ROOPSD mrv_histogram
#'
#' @examples
#' ## Generate sample
#' X = matrix( stats::rnorm( n = 100 ) , ncol = 4 )
#' Y = matrix( stats::rnorm( n = 100 ) , ncol = 4 )
#'
#' ## Fit dependence structure
#' ## Assume that the link beween column 2 and 4 is correct, and change also
#' ## the auto-correlation structure until lag 3 = lag_keep - 1
#' mvq = MVQuantilesShuffle$new( base::c(2,4) , lag_search = 6 , lag_keep = 4 )
#' mvq$fit(Y)
#' Z = mvq$transform(X)
#'
#' @export
MVQuantilesShuffle = R6::R6Class( "MVQuantilesShuffle" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field col_cond [vector] Conditionning columns
	col_cond   = NULL,
	#' @field col_ucond [vector] Un-conditionning columns
	col_ucond  = NULL,
	#' @field lag_search [integer] Number of lags to transform the dependence structure
	lag_search = NULL,
	#' @field lag_keep [integer] Number of lags to keep
	lag_keep   = NULL,
	
	#' @field n_features [integer] Number of features (dimensions), internal
	n_features = NULL,
	#' @field qY [matrix] Quantile structure fitted, internal
	qY         = NULL,
	#' @field bsYc [matrix] Block search fitted, internal
	bsYc       = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new MVQuantilesShuffle object.
    #' @param col_cond Conditionning colum
	#' @param lag_search  Number of lags to transform the dependence structure
	#' @param lag_keep Number of lags to keep
    #' @return A new `MVQuantilesShuffle` object.
	initialize = function( col_cond = base::c(1) , lag_search = 1 , lag_keep = 1 ) 
	{
		self$col_cond   = col_cond  
		self$lag_search = lag_search
		self$lag_keep   = lag_keep  
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit method
    #' @param Y [vector] Dataset to infer the dependance structure
    #' @return NULL
	fit = function(Y)
	{
		## Parameters
		self$n_features = ncol(Y)
		n_samplesY = nrow(Y)
		self$col_ucond = base::c()
		for( i in 1:self$n_features )
		{
			if( !(i %in% self$col_cond ) )
				self$col_ucond = base::c( self$col_ucond , i )
		}
		
		## Build non-parametric marginal distribution of Y
		rvY = mrv_histogram$new()$fit(Y)
		
		## Index to build block search matrix
		tiY = (n_samplesY + 1 - toeplitz(1:n_samplesY)[n_samplesY:1,])[1:self$lag_search,1:(n_samplesY-self$lag_search+1),drop=FALSE]
		
		## Find quantiles (i.e. ranks)
		self$qY  = rvY$cdf(Y)
		
		## Build conditionning block search
		qYc = self$qY[,self$col_cond,drop=FALSE]
		self$bsYc = base::c()
		for( i in 1:(n_samplesY-self$lag_search+1) )
		{
			self$bsYc = rbind( self$bsYc , as.vector(qYc[tiY[,i],]) )
		}
	},
	##}}}
	
	## transform ##{{{
	#' @description
    #' Transform method
    #' @param X [vector] Dataset to match the dependance structure with the Y fitted
    #' @return Z The X with the quantiles structure of Y
	transform = function(X)
	{
		## Parameters
		n_samplesX = nrow(X)
		
		## Build non-parametric marginal distribution of X
		rvX = mrv_histogram$new()$fit(X)
		
		## Index to build block search matrix
		tiX = (n_samplesX + 1 - toeplitz(1:n_samplesX)[n_samplesX:1,])[1:self$lag_search,1:(n_samplesX-self$lag_search+1),drop=FALSE]
		
		## Find quantiles (i.e. ranks)
		qX  = rvX$cdf(X)
		
		## Build conditionning block search
		## NOTE: in bsXc, the tiX[:,-1] column is added, otherwise the last values
		## are missing
		qXc = qX[,self$col_cond,drop=FALSE]
		bsXc = base::c()
		for( i in base::seq( 1 , n_samplesX-self$lag_search+1 , self$lag_keep ) )
		{
			bsXc = rbind( bsXc , as.vector(qXc[tiX[,i],]) )
		}
		bsXc = rbind( bsXc , as.vector(qXc[tiX[,ncol(tiX)],]) )
		
		## Now pairwise dist between cond. X / Y block search
		bsdistc = SBCK::pairwise_distances( bsXc , self$bsYc )
		idx_bsc = apply( bsdistc , 1 , which.min )
		
		## Find associated quantiles in unconditioning Y
		## NOTE: Here we split into lag_keep values, and some last missing values
		## lag_search - n_last is the numbers of last missing values.
		n_last = self$lag_search - (n_samplesX - (nrow(bsXc) - 1) * self$lag_keep) + 1
		
		## ===> Saved
#		qZuc = base::c()
#		for( i in idx_bsc[1:(length(idx_bsc)-1)] )
#		{
#			qZuc = rbind( qZuc , self$qY[,self$col_ucond,drop=FALSE][i:(i+self$lag_keep-1),,drop=FALSE] )
#		}
#		idxl = idx_bsc[length(idx_bsc)]
#		if( n_last < self$lag_search )
#			qZuc = rbind( qZuc , self$qY[,self$col_ucond,drop=FALSE][(idxl+n_last):(idxl+self$lag_search),,drop=FALSE] )
#		
#		## Now build qZ
#		qZ_unordered = cbind( qXc , qZuc )
#		qZ = matrix( 0 , nrow = nrow(qZ_unordered) , ncol = ncol(qZ_unordered) )
#		qZ[,base::c(self$col_cond,self$col_ucond)] = qZ_unordered
		## <===
		
		## ===> New
		qZ = base::c()
		for( i in idx_bsc[1:(length(idx_bsc)-1)] )
		{
			qZ = rbind( qZ , self$qY[i:(i+self$lag_keep-1),,drop=FALSE] )
		}
		idxl = idx_bsc[length(idx_bsc)]
		if( n_last < self$lag_search )
		{
			i0 = idxl + n_last
			i1 = idxl + self$lag_search
			if( i1 > nrow(self$qY) )
			{
				b  = i1 - nrow(self$qY)
				i0 = i0 - b
				i1 = i1 - b
			}
			qZ = rbind( qZ , self$qY[i0:i1,,drop=FALSE] )
		}
		qZ = qZ[1:nrow(X),]
		## <===
		
		## And finaly inverse quantiles
		Z = rvX$icdf(qZ)
		
		return(Z)
	}
	##}}}
	
	)
	
)
##}}}

## MVRanksShuffle ##{{{

#' MVRanksShuffle
#'
#' @description
#' Multivariate Schaake shuffle using the ranks.
#'
#' @details
#' Used to reproduce the dependence structure of a dataset to another dataset
#'
#' @references Vrac, M. et S. Thao (2020). “R2 D2 v2.0 : accounting for temporal
#'             dependences in multivariate bias correction via analogue rank
#'             resampling”. In : Geosci. Model Dev. 13.11, p. 5367-5387.
#'             doi :10.5194/gmd-13-5367-2020.
#'
#' @importFrom ROOPSD mrv_histogram
#'
#' @examples
#' ## Generate sample
#' X = matrix( stats::rnorm( n = 100 ) , ncol = 4 )
#' Y = matrix( stats::rnorm( n = 100 ) , ncol = 4 )
#'
#' ## Fit dependence structure
#' ## Assume that the link beween column 2 and 4 is correct, and change also
#' ## the auto-correlation structure until lag 3 = lag_keep - 1
#' mvr = MVRanksShuffle$new( base::c(2,4) , lag_search = 6 , lag_keep = 4 )
#' mvr$fit(Y)
#' Z = mvr$transform(X)
#'
#' @export
MVRanksShuffle = R6::R6Class( "MVRanksShuffle" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field col_cond [vector] Conditionning columns
	col_cond   = NULL,
	#' @field col_ucond [vector] Un-conditionning columns
	col_ucond  = NULL,
	#' @field lag_search [integer] Number of lags to transform the dependence structure
	lag_search = NULL,
	#' @field lag_keep [integer] Number of lags to keep
	lag_keep   = NULL,
	
	#' @field n_features [integer] Number of features (dimensions), internal
	n_features = NULL,
	#' @field qY [matrix] Ranks structure fitted, internal
	qY         = NULL,
	#' @field bsYc [matrix] Block search fitted, internal
	bsYc       = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new MVRanksShuffle object.
    #' @param col_cond Conditionning colum
	#' @param lag_search  Number of lags to transform the dependence structure
	#' @param lag_keep Number of lags to keep
    #' @return A new `MVRanksShuffle` object.
	initialize = function( col_cond = base::c(1) , lag_search = 1 , lag_keep = 1 ) 
	{
		self$col_cond   = col_cond  
		self$lag_search = lag_search
		self$lag_keep   = lag_keep  
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit method
    #' @param Y [vector] Dataset to infer the dependance structure
    #' @return NULL
	fit = function(Y)
	{
		## Parameters
		self$n_features = ncol(Y)
		n_samplesY = nrow(Y)
		self$col_ucond = base::c()
		for( i in 1:self$n_features )
		{
			if( !(i %in% self$col_cond ) )
				self$col_ucond = base::c( self$col_ucond , i )
		}
		
		## Index to build block search matrix
		tiY = (n_samplesY + 1 - toeplitz(1:n_samplesY)[n_samplesY:1,])[1:self$lag_search,1:(n_samplesY-self$lag_search+1),drop=FALSE]
		
		## Find quantiles (i.e. ranks)
		self$qY  = base::apply( Y , 2 , base::rank , ties.method = "first" )
		
		## Build conditionning block search
		qYc = self$qY[,self$col_cond,drop=FALSE]
		self$bsYc = base::c()
		for( i in 1:(n_samplesY-self$lag_search+1) )
		{
			self$bsYc = rbind( self$bsYc , as.vector(qYc[tiY[,i],]) )
		}
	},
	##}}}
	
	## transform ##{{{
	#' @description
    #' Transform method
    #' @param X [vector] Dataset to match the dependance structure with the Y fitted
    #' @return Z The X with the quantiles structure of Y
	transform = function(X)
	{
		## Parameters
		n_samplesX = nrow(X)
		n_dim      = ncol(X)
		
		## Index to build block search matrix
		tiX = (n_samplesX + 1 - toeplitz(1:n_samplesX)[n_samplesX:1,])[1:self$lag_search,1:(n_samplesX-self$lag_search+1),drop=FALSE]
		
		## Find quantiles (i.e. ranks)
		qX  = base::apply( X , 2 , base::rank , ties.method = "first" )
		
		## Shrink
		qY = self$qY
		n_samplesY = nrow(qY)
		if( n_samplesY < n_samplesX )
		{
			qX = round( qX * n_samplesY / n_samplesX )
		}
		if( n_samplesY > n_samplesX )
		{
			qY = round( qY * n_samplesX / n_samplesY )
		}
		
		## Build conditionning block search
		## NOTE: in bsXc, the tiX[:,-1] column is added, otherwise the last values
		## are missing
		qXc = qX[,self$col_cond,drop=FALSE]
		bsXc = base::c()
		for( i in base::seq( 1 , n_samplesX-self$lag_search+1 , self$lag_keep ) )
		{
			bsXc = rbind( bsXc , as.vector(qXc[tiX[,i],]) )
		}
		bsXc = rbind( bsXc , as.vector(qXc[tiX[,ncol(tiX)],]) )
		
		## Now pairwise dist between cond. X / Y block search
		bsdistc = SBCK::pairwise_distances( bsXc , self$bsYc )
		idx_bsc = apply( bsdistc , 1 , which.min )
		
		## Find associated quantiles in unconditioning Y
		## NOTE: Here we split into lag_keep values, and some last missing values
		## lag_search - n_last is the numbers of last missing values.
		n_last = self$lag_search - (n_samplesX - (nrow(bsXc)-1) * self$lag_keep)
		
		##===> Saved
#		qZuc = base::c()
#		for( i in idx_bsc[1:(length(idx_bsc)-1)] )
#		{
#			qZuc = rbind( qZuc , self$qY[,self$col_ucond,drop=FALSE][i:(i+self$lag_keep-1),,drop=FALSE] )
#		}
#		idxl = idx_bsc[length(idx_bsc)]
#		if( n_last < self$lag_search )
#			qZuc = rbind( qZuc , self$qY[,self$col_ucond,drop=FALSE][(idxl+n_last):(idxl+self$lag_search-1),,drop=FALSE] )
#		
#		## Now build qZ
#		qZ_unordered = cbind( qXc , qZuc )
#		qZ = matrix( 0 , nrow = nrow(qZ_unordered) , ncol = ncol(qZ_unordered) )
#		qZ[,base::c(self$col_cond,self$col_ucond)] = qZ_unordered
		##<===
		
		##===> New
		qZ = base::c()
		for( i in idx_bsc[1:(length(idx_bsc)-1)] )
		{
			qZ = rbind( qZ , self$qY[i:(i+self$lag_keep-1),,drop=FALSE] )
		}
		idxl = idx_bsc[length(idx_bsc)]
		if( n_last < self$lag_search )
			qZ = rbind( qZ , self$qY[(idxl+n_last):(idxl+self$lag_search-1),,drop=FALSE] )
		##<===
		
		## And finaly inverse ranks
		Xs = base::apply( X , 2 , sort )
		Z = NULL
		for( i in 1:n_dim )
		{
			Z = cbind( Z , Xs[qZ[,i],i] )
		}
		
		return(Z)
	}
	##}}}
	
	)
	
)
##}}}


