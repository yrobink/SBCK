
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

########################################################################
## WARNING                                                            ##
## Code here will be moved in future in an other package, called      ##
## ROOPSD: R Object Oriented Programming for Statistical Distribution ##
########################################################################

## rmultivariate_normal {{{

#' rmultivariate_normal
#'
#' Generate sample from a multivariate normal distribution
#'
#' @param n [integer]
#'        numbers of samples drawn
#' @param mean [vector]
#'        mean of Normal law
#' @param cov [matrix]
#'        covariance matrix
#'        
#' @return [matrix]
#'
#' @examples
#' X = SBCK::rmultivariate_normal( 10000 , base::c(0,0) , base::diag(2) )
#'
#' @export
rmultivariate_normal = function( n , mean , cov )
{
	svd = base::svd(cov)
	S =  svd$u %*% base::diag(base::sqrt(svd$d)) %*% base::t(svd$u)
	n_dim = length(mean)
	X = base::matrix( stats::rnorm(n*n_dim) , nrow = n , ncol = n_dim )
	X = base::t(base::apply( X , 1 , function(x) { mean + S %*% x } ))
#	X = base::matrix( NA , nrow = n , ncol = n_dim )
#	for( i in 1:n )
#	{
#		X[i,] = mean + S %*% stats::rnorm( n_dim , mean = 0 , sd = 1 )
#	}
	return(X)
}
##}}}

## rorthogonal_group {{{

#' rorthogonal_groupd
#'
#' Generate sample from a the orthogonal group O(d)
#'
#' @param d [integer]
#'        Dimension of the matrix
#' @param n [integer]
#'        numbers of samples drawn
#'        
#' @return [array], dim = d*d*n
#'
#' @examples
#' M = SBCK::rorthogonal_group( 2 , 10 )
#'
#' @export
rorthogonal_group = function( d , n = 1 )
{
	rot = array( stats::runif( d * d * n ) , dim = base::c( d , d , n ) )
	for( i in 1:n )
	{
		QR = base::qr( rot[,,i] )
		Q  = base::qr.Q(QR)
		R  = base::diag(base::qr.R(QR))
		R  = base::diag( R / base::abs(R) )
		rot[,,i] = Q %*% R
	}
	return(rot)
}
##}}}

## Abstract dist ##{{{

#' AbstractDist 
#'
#' Base class for OOP statistical distribution
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param ddist [function]
#'        Density function, e.g. dnorm
#' @param pdist [function]
#'        Distribution function, e.g. pnorm
#' @param qdist [function]
#'        Quantile function, e.g. qnorm
#' @param rdist [function]
#'        Random generator function, e.g. rnorm
#' @param freeze [boolean]
#'        If we freeze or note the distribution. TRUE after fit in general
#' @param n  [integer]
#'        For rvs function, number of samples drawn
#' @param x  [integer]
#'        For density and logdensity function. Vector of quantiles
#' @param q  [integer]
#'        For cdf and sf functions. Vector of quantiles
#' @param p  [integer]
#'        For icdf and isf functions. Vector of probabilities
#' @param Y  [integer]
#'        For fit function. Dataset to infer parameters. (max likelihood)
#'
#' @return Object of \code{\link{R6Class}} with methods to use a statistical distribution.
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(ddist,pdist,qdist,rdist,freeze)}}{This method is used to create object of this class with \code{AbstractDist}}
#'   \item{\code{rvs(n)}}{Draw n samples}.
#'   \item{\code{density(x)}}{Density along vector of quantile x}.
#'   \item{\code{logdensity(x)}}{log of density along vector of quantile x}.
#'   \item{\code{cdf(q)}}{Cumulative Distribution Function along vector of quantile q}.
#'   \item{\code{sf(q)}}{Survival function (1-CDF) along vector of quantile q}.
#'   \item{\code{icdf(p)}}{Inverse of cdf along vector of probabilities p}.
#'   \item{\code{isf(p)}}{Inverse of sf along vector of probabilities p}.
#'   \item{\code{fit(Y)}}{Fit function to infer parameters. By max likelihood}.
#' }
#' @examples
#' ##
#' ##
#' @export
AbstractDist = R6::R6Class( "AbstractDist",
	
	## Public elements
	##============={{{
	public = list(
	
	## Arguments
	##==========
	ddist = NULL,
	pdist = NULL,
	qdist = NULL,
	rdist = NULL,
	freeze = FALSE,
	
	## Constructor
	##============
	
	initialize = function( ddist , pdist , qdist , rdist , freeze = FALSE , parametric = TRUE )##{{{
	{
		print("=====================\nNot call!!!!!\n=====================\n")
		self$ddist = ddist
		self$pdist = pdist
		self$qdist = qdist
		self$rdist = rdist
		self$freeze = freeze
		private$.is_parametric = parametric
	},
	##}}}
	
	## Methods
	##========
	
	rvs = function( n )##{{{
	{
		params = private$params()
		params$n = n
		return(base::do.call( self$rdist , params ))
	},
	##}}}
	
	density = function(x)##{{{
	{
		params = private$params()
		params$x = x
		return(base::do.call( self$ddist , params ) )
	},
	##}}}
	
	logdensity = function(x)##{{{
	{
		params = private$params()
		params$x = x
		params$log = TRUE
		return(base::do.call( self$ddist , params ) )
	},
	##}}}
	
	cdf = function(q)##{{{
	{
		params = private$params()
		params$q = q
		return(base::do.call( self$pdist , params ) )
	},
	##}}}
	
	sf = function(q)##{{{
	{
		params = private$params()
		params$q = q
		params$lower.tail = FALSE
		return(base::do.call( self$pdist , params ) )
	},
	##}}}
	
	icdf = function(p)##{{{
	{
		params = private$params()
		params$p = p
		params$lower.tail = TRUE
		return(base::do.call( self$qdist , params ) )
	},
	##}}}
	
	isf = function(p)##{{{
	{
		params = private$params()
		params$p = p
		params$lower.tail = FALSE
		return(base::do.call( self$qdist , params ) )
	},
	##}}}
	
	is_parametric = function()##{{{
	{
		return(private$.is_parametric)
	},
	##}}}
	
	## Fit
	##====
	
	fit = function(Y)##{{{
	{
		private$fit_initialization(Y)
		opt = stats::optim( par = as.vector(private$params()) , fn = private$negloglikelihood , method = "BFGS" , Y = Y )
		private$set_params(opt$par)
	}
	##}}}
	
	),
	##}}}
	
	## Private elements
	##==============={{{
	
	private = list(
	
	.is_parametric = NULL,
	
	params = function()##{{{
	{},
	##}}}
	
	set_params = function(params)##{{{
	{},
	##}}}
	
	negloglikelihood = function( params , Y )##{{{
	{
		private$set_params(params)
		return(-base::sum(self$logdensity(Y)))
	},
	##}}}
	
	fit_initialization = function(Y)##{{{
	{}
	##}}}
	
	)
	##}}}
	
)
##}}}

## ROOPSD_Normal ##{{{

#' ROOPSD_Normal 
#'
#' Normal distribution in OOP way. Based on AbstractDist
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param mean [function]
#'        Mean
#' @param sd [function]
#'        Standard deviation
#' @param freeze [boolean]
#'        If we freeze or note the distribution. TRUE after fit in general
#'
#' @return Object of \code{\link{R6Class}} with methods to use a statistical distribution.
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(mean,sd,freeze)}}{This method is used to create object of this class with \code{ROOPSD_Normal}}
#' }
#' @examples
#' ##
#' ##
#' @export
ROOPSD_Normal = R6::R6Class( "ROOPSD_Normal",
	
	inherit = AbstractDist,
	
	## Public elements
	##============={{{
	
	public = list(
	
	## Arguments
	##==========
	mean = 0.,
	sd   = 1.,
	
	
	## Constructor
	##============
	
	initialize = function( mean = 0 , sd = 1 , freeze = FALSE )##{{{
	{
		super$initialize( stats::dnorm , stats::pnorm , stats::qnorm , stats::rnorm , freeze )
		self$mean = mean
		self$sd   = sd
	}
	##}}}
	
	),
	##}}}
	
	## Private elements
	##==============={{{
	private = list(
	
	params = function()##{{{
	{
		return( list( mean = self$mean , sd = self$sd ) )
	},
	##}}}
	
	set_params = function(params)##{{{
	{
		self$mean = params[1]
		self$sd   = params[2]
	},
	##}}}
	
	fit_initialization = function(Y)##{{{
	{
		self$mean = base::mean(Y)
		self$sd   = stats::sd(Y)
	}
	##}}}
	
	)
	##}}}
)
##}}}

## ROOPSD_Exponential ##{{{

#' ROOPSD_Exponential 
#'
#' Exponential distribution in OOP way. Based on AbstractDist
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param Rate [function]
#'        Rate = 1 / mean
#' @param freeze [boolean]
#'        If we freeze or note the distribution. TRUE after fit in general
#'
#' @return Object of \code{\link{R6Class}} with methods to use a statistical distribution.
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(rate,freeze)}}{This method is used to create object of this class with \code{ROOPSD_Exponential}}
#' }
#' @examples
#' ##
#' ##
#' @export
ROOPSD_Exponential = R6::R6Class( "ROOPSD_Exponential",
	
	inherit = AbstractDist,
	
	## Public elements
	##============={{{
	
	public = list(
	
	## Arguments
	##==========
	rate = 1.,
	
	
	## Constructor
	##============
	
	initialize = function( rate = 1. , freeze = FALSE )##{{{
	{
		super$initialize( stats::dexp , stats::pexp , stats::qexp , stats::rexp , freeze )
		self$rate = rate
	}
	##}}}
	
	),
	##}}}
	
	## Private elements
	##==============={{{
	private = list(
	
	params = function()##{{{
	{
		return( list( rate = self$rate ) )
	},
	##}}}
	
	set_params = function(params)##{{{
	{
		self$rate = params[1]
	},
	##}}}
	
	fit_initialization = function(Y)##{{{
	{
		self$rate = 1. / base::mean(Y)
	}
	##}}}
	
	)
	##}}}
)
##}}}

## ROOPSD_Gamma ##{{{

#' ROOPSD_Gamma 
#'
#' Gamma distribution in OOP way. Based on AbstractDist
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param shape [function]
#'        Shape
#' @param scale [function]
#'        scale
#' @param freeze [boolean]
#'        If we freeze or note the distribution. TRUE after fit in general
#'
#' @return Object of \code{\link{R6Class}} with methods to use a statistical distribution.
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(rate,freeze)}}{This method is used to create object of this class with \code{ROOPSD_Gamma}}
#' }
#' @examples
#' ##
#' ##
#' @export
ROOPSD_Gamma = R6::R6Class( "ROOPSD_Gamma",
	
	inherit = AbstractDist,
	
	## Public elements
	##============={{{
	
	public = list(
	
	## Arguments
	##==========
	shape = 0.,
	scale = 1.,
	
	
	## Constructor
	##============
	
	initialize = function( shape = 0.5 , scale = 1 , freeze = FALSE )##{{{
	{
		super$initialize( stats::dgamma , stats::pgamma , stats::qgamma , stats::rgamma , freeze )
		self$shape = shape
		self$scale = scale
	}
	##}}}
	
	),
	##}}}
	
	## Private elements
	##==============={{{
	private = list(
	
	params = function()##{{{
	{
		return( list( shape = self$shape , scale = self$scale ) )
	},
	##}}}
	
	set_params = function(params)##{{{
	{
		self$shape = params[1]
		self$scale = params[2]
	},
	##}}}
	
	fit_initialization = function(Y)##{{{
	{
		e = base::mean(Y)
		v = stats::var(Y)
		self$shape = e^2 / v
		self$scale = v / e
	}
	##}}}
	
	)
	##}}}
)
##}}}

## ROOPSD_rv_histogram ##{{{

#' ROOPSD_rv_histogram 
#'
#' Empirical distribution in OOP way. Use quantile.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param n  [integer]
#'        For rvs function, number of samples drawn
#' @param x  [integer]
#'        For density and logdensity function. Vector of quantiles
#' @param q  [integer]
#'        For cdf and sf functions. Vector of quantiles
#' @param p  [integer]
#'        For icdf and isf functions. Vector of probabilities
#' @param Y  [vector]
#'        For fit function. Dataset to infer parameters. (Quantile)
#' @param bins [integer]
#'        Numbers of quantiles to estimate between 0 and 1.
#' @return Object of \code{\link{R6Class}} with methods to use a statistical distribution.
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{ROOPSD_rv_histogram}}
#'   \item{\code{rvs(n)}}{Draw n samples}.
#'   \item{\code{cdf(q)}}{Cumulative Distribution Function along vector of quantile q}.
#'   \item{\code{sf(q)}}{Survival function (1-CDF) along vector of quantile q}.
#'   \item{\code{icdf(p)}}{Inverse of cdf along vector of probabilities p}.
#'   \item{\code{isf(p)}}{Inverse of sf along vector of probabilities p}.
#'   \item{\code{fit(Y,bins)}}{Fit function to infer parameters}.
#' }
#' @examples
#' ##
#' ##
#' @export
ROOPSD_rv_histogram = R6::R6Class( "ROOPSD_rv_histogram" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	min   = NULL,
	max   = NULL,
	freeze = FALSE,
	tol = 1e-2,
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{},
	
	rvs = function( n )
	{
		p = stats::runif( n , 0 , 1 )
		return(self$icdf(p))
	},
	
	cdf = function( q )
	{
		return(private$cdffn(q))
	},
	
	icdf = function( p )
	{
		return(private$icdffn(p))
	},
	
	sf = function( q )
	{
		return( 1. - private$cdffn(q) )
	},
	
	isf = function( p )
	{
		return(private$icdffn(1. - p))
	},
	
	is_parametric = function()
	{
		return(FALSE)
	},
	
	fit = function( Y , bins = as.integer(100) )
	{
		self$min = base::min(Y)
		self$max = base::max(Y)
		
		## CDF and iCDF function
		private$cdffn = stats::ecdf(Y)
		x = NULL
		delta = NULL
		if( is.integer(bins) )
		{
			delta = 1e-2 * (self$max - self$min)
			x = base::seq( self$min - delta , self$max + delta , length = bins )
		}
		else
		{
			x = bins
			delta = min( diff(bins) )
		}
		quants = private$cdffn(x)
		private$icdffn = stats::approxfun( quants , x , yleft = self$min - delta , yright = self$max + delta )
	}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	cdffn  = NULL,
	icdffn = NULL
	
	
	)
)
##}}}

## ROOPSD_rv_ratio_histogram ##{{{

#' ROOPSD_rv_ratio_histogram
#'
#' Empirical distribution for ratio variables in OOP way. We fit separatly P( X < x | X > 0 ) and P(X=0)
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param n  [integer]
#'        For rvs function, number of samples drawn
#' @param x  [integer]
#'        For density and logdensity function. Vector of quantiles
#' @param q  [integer]
#'        For cdf and sf functions. Vector of quantiles
#' @param p  [integer]
#'        For icdf and isf functions. Vector of probabilities
#' @param Y  [vector]
#'        For fit function. Dataset to infer parameters. (Quantile)
#' @param bins [integer]
#'        Numbers of quantiles to estimate between 0 and 1.
#' @return Object of \code{\link{R6Class}} with methods to use a statistical distribution.
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{ROOPSD_rv_ratio_histogram}}
#'   \item{\code{rvs(n)}}{Draw n samples}.
#'   \item{\code{cdf(q)}}{Cumulative Distribution Function along vector of quantile q}.
#'   \item{\code{sf(q)}}{Survival function (1-CDF) along vector of quantile q}.
#'   \item{\code{icdf(p)}}{Inverse of cdf along vector of probabilities p}.
#'   \item{\code{isf(p)}}{Inverse of sf along vector of probabilities p}.
#'   \item{\code{fit(Y,bins)}}{Fit function to infer parameters}.
#' }
#' @examples
#' ##
#' ##
#' @export
ROOPSD_rv_ratio_histogram = R6::R6Class( "ROOPSD_rv_ratio_histogram" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	distP = NULL,
	p0    = 0,
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{},
	
	rvs = function( n )
	{
		p = stats::runif( n , 0 , 1 )
		return(self$icdf(p))
	},
	
	cdf = function( q )
	{
		cdf = base::rep( 0 , length(q) )
		idxp = q > 0
		idx0 = !idxp
		cdf[idxp] = ( 1 - self$p0 ) * self$distP$cdf(q[idxp]) + self$p0
		cdf[idx0] = self$p0 / 2
		return(cdf)
	},
	
	icdf = function( p )
	{
		idxp = p > self$p0
		idx0 = !idxp
		icdf = base::rep( 0 , length(p) )
		icdf[idxp] = self$distP$icdf( (p[idxp] - self$p0) / ( 1 - self$p0 ) )
		icdf[idx0] = 0
		return(icdf)
	},
	
	sf = function( q )
	{
		return( 1. - self$cdf(q) )
	},
	
	isf = function( p )
	{
		return(self$icdf(1. - p))
	},
	
	is_parametric = function()
	{
		return(FALSE)
	},
	
	fit = function( X , bins = 100 )
	{
		Xp = X[X>0]
		self$distP = rv_ratio_histogram$new()
		self$distP$fit( Xp , bins = bins )
		self$p0 = base::sum(!(X>0)) / length(X)
	}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	cdffn  = NULL,
	icdffn = NULL
	
	
	)
)
##}}}

