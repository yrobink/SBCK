
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

