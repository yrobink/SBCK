
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

#' Shift
#'
#' Class to transform a dataset such that autocorrelations are intervariables correlations
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param lag [integer]
#'        max lag for autocorrelations
#' @param method [character]
#'        If "row" inverse by row, else by column
#' @param ref [integer]
#'        starting point for inverse transform
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(lag,method,ref,)}}{This method is used to create object of this class with \code{Shift}}
#'   \item{\code{transform(X)}}{Method to shift a dataset}
#'   \item{\code{inverse(Xs)}}{Method to inverse the shift of a dataset}
#' }
#' @examples
#' X = base::t(matrix( 1:20 , nrow = 2 , ncol = 10 ))
#' 
#' sh = Shift$new(1)
#' Xs = sh$transform(X)
#' Xi = sh$inverse(Xs)
#'
#' @export
Shift = R6::R6Class( "Shift" ,
	
	private = list( ##{{{
	
	.method = NULL,
	.ref    = NULL,
	
	inverse_by_row = function(Xs)##{{{
	{
		n_features = as.integer( dim(Xs)[2] / ( self$lag + 1 ))
		n_samples  = dim(Xs)[1] + self$lag
		Xi = matrix( NA , ncol = n_features , nrow = n_samples )
		for( r in base::c(1:self$lag,self$ref) )
		{
			idx = base::seq( r , n_samples - self$lag , self$lag )
			Xs0 = as.vector(base::t(Xs[idx[-length(idx)],1:(base::ncol(Xs) - n_features)]))
			Xs1 = as.vector(base::t(Xs[idx[length(idx)],]))
			Xs01 = base::t(matrix( base::c(Xs0,Xs1) , nrow = n_features ))
			Xi[r:(r+base::nrow(Xs01)-1),] = Xs01
		}
		
		return(Xi)
	},
	##}}}
	
	inverse_by_col = function(Xs)##{{{
	{
		n_features = as.integer( dim(Xs)[2] / (self$lag + 1) )
		n_samples  = dim(Xs)[1] + self$lag
		Xi = matrix( NA , nrow = n_samples , ncol = n_features )
		
		ref = self$ref
		if( ref < 1 )
			ref = 1
		if( ref > n_features )
			ref = n_features
		
		for( i in base::c( 0:self$lag , ref ) )
		{
			db = i * n_features + 1
			de = i * n_features + n_features
			tb = i + 1
			te = n_samples - ( self$lag + 1 ) + i + 1
			Xi[tb:te,] = Xs[,db:de]
		}
		
		return(Xi)
	}
	##}}}
	
	),
	##}}}
	
	active = list( ##{{{
	
	method = function(value)##{{{
	{
		if( missing(value) )
		{
			return(private$.method)
		}
		else
		{
			if( value %in% base::c("col","row") )
				private$.method = value
			else
				private$.method = "row"
		}
	},
	##}}}
	
	ref = function(value)##{{{
	{
		if( missing(value) )
		{
			return(private$.ref)
		}
		else
		{
			private$.ref = ( (value - 1) %% ( self$lag + 1 ) ) + 1
		}
	}
	
	),
	##}}}
	
	##}}}
	
	public = list( ##{{{
	
	###############
	## Arguments ##
	###############
	
	lag    = NULL,
	
	#################
	## Constructor ##
	#################
	
	initialize = function( lag , method = "row" , ref = 1 )##{{{
	{
		self$lag    = lag
		self$ref    = ref
		self$method = method
	},
	##}}}
	
	
	#############
	## Methods ##
	#############
	
	transform = function( X )##{{{
	{
		if( is.vector(X) )
			X = matrix( X , nrow = length(X) , ncol = 1 )
		
		n_samples  = base::dim(X)[1]
		n_features = base::dim(X)[2]
		Xs = matrix( NA , nrow = n_samples - self$lag , ncol = ( self$lag + 1 ) * n_features )
		
		for( i in 0:self$lag )
		{
			db = i * n_features + 1
			de = i * n_features + n_features
			tb = i + 1
			te = n_samples - ( self$lag + 1 ) + i + 1
			Xs[,db:de] = X[tb:te,]
		}
		
		return(Xs)
	},
	##}}}
	
	inverse = function( Xs )##{{{
	{
		if( self$method == "row" )
			return(private$inverse_by_row(Xs))
		else
			return(private$inverse_by_col(Xs))
	}
	##}}}
	
	)
	##}}}
)

