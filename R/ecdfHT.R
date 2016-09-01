#' The \code{ecdfHT} package computes and plot a transformed empirical cdf for data.
#' This is useful because a standard empirical cdf (ecdf) gives little information
#' about the tails of the data when there are extreme values.
#'
#' The transform is nonparametric: linear in the middle of the data and matched to a log-log
#' transform on the tails, where the tail regions  are determined by quantiles.
#' If the data has power law behavior on the tails, the plot is linear on those tails,
#' so this plot can be used as a graphical diagnostic to determine if a data set is heavy tailed.
#'
#' In addition, there are functions to
#' \itemize{
#' \item annotate the plot, add custom axes and grid lines
#' \item overlay proposed models on the plot
#' \item fit the tails using linear regression on the transformed tails
#' \item combine the empirical cdf in the middle and the above fit on the tails to get a semi-parametric fit to the data
#' \item compute cdf, pdf, quantiles, and simulate from the semi-parametric fit
#' \item some multivariate plots that look at tail behavior of multiple components and some idea of the dependence.
#' }
#'
#' I will try to fix the code if you provide a simple demonstration of a bug.  Polite suggestions for
#' improvements will be considered if there is time available.
#'
#' @seealso  \code{\link{ecdfHT}} for a basic plot,
#' \code{\link{ecdfHT.draw}} for annotations and additions to a basic plot,
#' \code{\link{ecdfHT.fit}} to fit a semi-parametric model to the data,
#' \code{\link{pecdfHT}} to compute the cdf, pdf, quantiles and simulate from a
#' semi-parametric model,
#' \code{\link{ecdfHT.multivar}} for multivariate generalizations.
#'
#'
#' @importFrom "grDevices"  "dev.new"
#' @importFrom "graphics" "abline"
#' @importFrom "graphics" "legend"
#' @importFrom "graphics" "lines"
#' @importFrom "graphics" "mtext"
#' @importFrom "graphics"  "par"
#' @importFrom "graphics" "plot"
#' @importFrom "graphics"  "points"
#' @importFrom "graphics" "title"
#' @importFrom "graphics" "pairs"
#' @importFrom "stats" "IQR"
#' @importFrom "stats" "lm"
#' @importFrom "stats" "median"
#' @importFrom "stats" "quantile"
#' @importFrom "stats" "runif"
#' @importFrom "stats" "approxfun"
#'
#' @title ecdfHT: A package to plot an empirical cdf for heavy tailed data
#' @docType package
#' @name ecdfHT-package
NULL
###################################################################################
#' Plot a transformed empirical cdf
#'
#' Produces a basic plot showing a transformed empirical cdf for heavy tailed data.  It uses a log-log
#' transform on the tails, which shows power law decay as linear.
#'
#' @param x A vector of data
#' @param scale.q A vector of 3 probabilites; specifies the quantiles of the data to use for the left tail, mid region, and right tail
#' @param show.axes.labels Boolean value: indicates whether default labels are plotted or not.  (Use function \code{ecdfHT.axes} to add custom labels)
#' @param show.plot Boolean value: show plot or only do calculations
#' @param type Type of plot, passed to plot.  Use type='p' for points, type='l' for lines
#' @param ... Optional graphical parameters, e.g. col='red'
#'
#' @details Most of the work is done by \code{ecdfHT.draw} and the associated helper functions.
#'
#' Assuming no repeats in x,  ecdf = (standard ecdf - (1/2))/n, like type=5 in the R function \code{quantile}.
#' So instead of taking values 1/n, 2/n, 3/n, ... , k/n, ...,  1 it takes values
#'   1/(2n), 3/(2n), ..., (2k-1)/(2n), ..., (2n-1)/(2n).
#' This avoids 0 at lower endpoint and 1 at upper endpoint, which causes problems
#' when we extend tails with a power law.  (If there are m repeated x values, then the corresponding
#' jump in the ecdf at that point is m/n instead of 1/n.)
#'
#' @return An object of class 'ecdfHT.transform' which gives the information necessary to draw the plot and later add other curves and labels.
#' This list is returned invisibly and contains the following fields:
#' \describe{
#' \item{scale.q}{vector of length 3, copied from the input argument}
#' \item{scale.x}{vector of length 3, the quantiles from the data corresponding to scale.q}
#' \item{xsort}{vector of the sorted, unique data values}
#' \item{ecdf}{nonstandard empirical cdf, see details}
#' \item{xx}{transformed x values: xx[i]=h(xsort[i])}
#' \item{yy}{transformed ecdf values: yy[i]=g(ecdf[i])}
#' }
#'
#' @export
#' @details The default values scale.q=c(.25,.5,.75) splits the data into quartiles; picking
#' different quantiles splits the data into 4 different groups: the lowest group is the left tail,
#' i.e. all values less than the quantile corresponding to scale.q[1];
#' the next group is between the left tail and center = quantile scale.q[2]); the third
#' group is the center and quantile scale.q[3]; the last group is the upper tail.
#' For two-sided data, it makes sense to use something like (p,0.5,1-p) for scale.q, where
#' p is choosen to determine where the tail regions begin.
#'
#' For one-sided data, it makes sense to use scale.q=c(0,0,p).  In this case,
#' the first two groups are empty and the effect is to divide the data into two groups:
#' a moderate/lower range and a right tail.  See the example below with nonnegative data.
#'
#' The transformations h(x) acts on these different regions.  It is linear on the middle
#' two regions and logarithmic on the tails.  The transformation g(p) acts on the corresponding
#' values of the ecdf described above.
#' The basic plot shows (h(x[i]),g(ecdf[i])): the first component is a monotonic transform of the
#' x values, the second component is a monotonic transform of the ecdf.
#' See the accompanying vignette for exact definitions: go to the package index and click on User
#' guides, package vignettes and other documentation.
#'
#' @seealso  \code{\link{ecdfHT.draw}} for annotations and additions to a basic plot
#'
#' @examples
#' x <- rcauchy( 1000 )
#' ecdfHT( x )
#' title("basic ecdfHT plot")
#'
#' xabs <- abs(x)
#' ecdfHT( xabs, scale.q=c(0,0,.75) )
#' title("one sided data")
#'
ecdfHT <- function( x, scale.q=c(0.25,0.5,0.75), show.axes.labels=TRUE, show.plot=TRUE, type='p', ... ) {

stopifnot( is.vector(x), is.vector(scale.q), is.numeric(x), is.numeric(scale.q),
           is.logical(show.axes.labels), is.logical(show.plot),
           length(scale.q)==3, scale.q[1] <= scale.q[2], scale.q[2] <= scale.q[3] )

# compute cut points for the transform of the variables
scale.x <- quantile(x,prob=scale.q,names=FALSE)

# save transform info for later use
transform.info <- list( scale.q=scale.q, scale.x=scale.x)

# sort x and select unique values; then compute "non-standard" ecdf, handling duplicates.
xsort <- unique(sort(x))
ecdf <- ( cumsum(tabulate(match(x, xsort))) - 0.5 )/length(x)

r <- ecdfHT.draw(transform.info,xsort,ecdf,show.plot=show.plot,new.plot=TRUE,show.ci=FALSE,type=type,...)

if(show.plot) {
  # draw main plot
  if (show.axes.labels) {
    # default labels and style for the graph axes
    y.labels <- c(.01,.05,.1,.5,.90,.95,.99)
    x.labels <- quantile(x, c(.01,.1,.5,.9,.99) )
    x.labels <- signif( x.labels, 3 )  # round to show 3 sig. digits
    ecdfHT.axes( transform.info, x.labels, y.labels, FALSE, FALSE )
  }
}

# add more fields to the transform info
transform.info$xsort <- xsort
transform.info$ecdf <- ecdf
transform.info$xx <- r$xx
transform.info$yy <- r$yy
class(transform.info) <- "ecdfHT.transform"

invisible(transform.info) }
#######################################################################
#' Graph and annotate an ecdfHT plot
#'
#' Does the computations and plotting for \code{ecdfHT} and can be used to add to an existing plot.
#' @param transform.info A list with information about the transformation, computed in \code{ecdfHT}
#' @param x The data, a vector of double precision numbers.  Assumbed to be sorted and have distinct values.
#' @param p Probabilities, a vector of doubles.  Typically p[i]=(i=0.5)/length(x), unless there are repeats in x.
#' @param show.plot Boolean value: indicates whether to plot or not.
#' @param new.plot Boolean value: indicates whether to produce a new plot or add to an existing plot.
#' @param show.ci Boolean value: indicates whether or not confidence intervals are shown.
#' @param xlab String to label the horizontal axis.
#' @param ylab String to label the vertical axis.
#' @param ... Optional parameters for the plot, e.g. col='red'.
#'
#' @details  \code{ecdfHT.draw} computes transform and plots.
#' \code{ecdfHT.axes} draws axes on the plot; it can be used to manually select tick marks, etc.
#' \code{ecdfHT.h} computes the function h(x) for the transformation of the horizontal axis.
#' \code{ecdfHT.g} computes the function g(p) for the transformation of the vertical axis.
#'
#' Always call \code{ecdfHT} first to produce the basic plot, then use \code{ecdfHT.draw}
#' to add other curves to the plot as in the examples below
#'
#' @return A list of values used in the plot, see return value of \code{ecdfHT}.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- rcauchy( 1000 )
#' t.info <- ecdfHT( x, show.axes=FALSE )
#' ecdfHT.axes( t.info, x.labels=c(-50,-5,0,5,50), y.labels=c(.001,.01,.1,.5,.9,.99,.999),
#'   show.vert.gridlines=TRUE, show.horiz.gridline=TRUE, lty=2 )
#' q1 <- qcauchy(t.info$ecdf) # Cauchy quantiles
#' ecdfHT.draw( t.info, q1, t.info$ecdf, col='red',show.ci=TRUE)
#' q2 <- qnorm(t.info$ecdf,sd=sd(x))  # Gaussian quantiles
#' ecdfHT.draw( t.info, q2, t.info$ecdf, col='green',show.ci=TRUE)
#' title(paste("simulated Cauchy data, n=",length(x),"\nred=Cauchy cdf, green=normal cdf"))
#'
ecdfHT.draw <- function( transform.info, x, p, show.plot=TRUE, new.plot=FALSE, show.ci=FALSE, xlab="x", ylab="", ...) {

# transform x and p values
xx <- ecdfHT.h( x, transform.info$scale.x )
yy <- ecdfHT.g( p, transform.info$scale.q )

# draw plot
if (show.plot) {
  if (new.plot) {  # do basic plot
    par(xaxt="n",yaxt="n")
    plot(xx,yy, xlab=xlab, ylab=ylab, ...)
  } else {  # add to existing plot
    lines(xx,yy,...)
  }

  # show pointwise confidence intervals if requested
  if (show.ci) {
    n <- length(x)
    low <- pmax(1e-100,p - 1.96*sqrt(p*(1-p)/n) )
    low <- ecdfHT.g(low,transform.info$scale.q)
    lines(xx,low,lty=3,...)
    hi <- pmin(1,p + 1.96*sqrt(p*(1-p)/n))
    hi <- ecdfHT.g(hi, transform.info$scale.q)
    lines(xx,hi,lty=3,...)
  }
}
invisible(list(xx=xx,yy=yy,x=x,p=p)) }
#######################################################################
#' @rdname ecdfHT.draw
#' @export
#' @param x.labels Vector of numbers specifying the location of the labels on the horizontal axis
#' @param y.labels Vector of numbers specifying the location of the labels on the vertical axis
#' @param show.vert.gridlines Boolean value indicating whether or not vertical grid lines should be drawn.
#' @param show.horiz.gridlines Boolean value indicating whether or not horizontal grid lines should be drawn.
#'
#'
ecdfHT.axes <- function( transform.info, x.labels=c(), y.labels= c(),
                         show.vert.gridlines=FALSE, show.horiz.gridlines=FALSE, ... ) {
# label axes with nonlinear scaling
# use ... arguments to control line type and color of gridlines

# could not get axis( ) command to work as desired, draw things ourselves
par.usr <- par("usr") # left, right, bottom, top user coordinates
tick.width <- 0.01*(par.usr[2]-par.usr[1])
tick.height <- 0.01*(par.usr[4]-par.usr[3])
x.loc <- NULL; y.loc <- NULL

if( length(x.labels)>0 ){
  # compute x labels and plot
  x.loc <- ecdfHT.h(x.labels,transform.info$scale.x)
  x.str <- paste( signif(x.labels,4))
  mtext( x.str,side=1,at=x.loc,adj=0.5)
  for (i in 1:length(x.labels) ) {
    if (show.vert.gridlines) {
      abline(v=x.loc[i],...)  # vertical grid lines
    } else {
      lines(rep(x.loc[i],2),c(par.usr[3],par.usr[3]+tick.height),xpd=NA ) # tick marks only
    }
  }
}
if( length(y.labels)>0 ){
  # compute y labels and plot
  y.loc <- ecdfHT.g(y.labels,transform.info$scale.q )
  y.str <- paste(y.labels)
  for (i in 1:length(y.str)) {  # avoid exponential format for moderately small probabilities
    if ((y.labels[i] < .001) && (y.labels[i] > .000001)) {
      y.str[i] <- formatC(y.labels[i],digits=ceiling(-log(y.labels[i],base=10)))
    }
  }
  mtext(paste(y.str," "),side=2,at=y.loc,padj=0.5,las=1)
  for (i in 1:length(y.labels) ) {
    if (show.horiz.gridlines) {
      abline(h=y.loc[i],...)  # horizontal grid lines
    } else {
      lines( par.usr[1]+tick.width*c(0,1), rep(y.loc[i],2),xpd=NA )  # tick marks only
    }
  }
}
invisible( list(x.labels=x.labels,x.loc=x.loc,y.labels=y.labels,y.loc=y.loc))}
####################################################################
#' @rdname ecdfHT.draw
#' @param t A vector of length 3 that specifies the x values that determine the left tail, middle, and right tail
#' @export
#' @return \code{ecdfHT.h} returns the vector y=h(x;t), \code{ecdfHT.g} returns the vector y=g(p;q)
#' @examples
#' x <- seq(-5,5,1)
#' t <- c(-3,0,3)
#' ecdfHT.h(x,t)
#' p <- seq(0.05,.95,.1)
#' q <- c(.1,.5,.9)
#' ecdfHT.g(p,q)
#'
ecdfHT.h <- function(x,t) {
# compute the transformation h(x) on the scale of the data

# shift and scale the data
y <- x
i <- which(x >= t[2])
y[i] <- (x[i]-t[2])/(t[3]-t[2])
i <- which(x < t[2])
y[i] <- (x[i]-t[2])/(t[2]-t[1])

# now compute log transform on lower and then upper tail
lo <- which(y < -1)
y[lo] <- -1 - log( -y[lo] )

hi <- which(y > 1)
y[hi] <- 1 + log( y[hi] )

return(y) }
#######################################################################
#' @rdname ecdfHT.draw
#' @param q A vector of length 3 that specifies the quantile values that determine the left tail, middle, and right tail.
#' @export
#'
ecdfHT.g <- function(p,q) {
# compute the transformation g(x) on the probability scale (0,1)

y <- p
lo <- which(p < q[1])
y[lo] <- q[1]*(1 + log( p[lo]/q[1] ) )

hi <- which(p > q[3])
y[hi] <- q[3] + (q[3]-1)* log( (1.0 -p[hi])/(1-q[3]) )

return(y) }
########################################################################
#' Fit heavy tailed data with a semi-parameteric model
#'
#' Compute an interpolation of the transformed cdf in the middle with parametric power law decay on the tails.
#' @param p Vector of 2 probabilities that identify the quantile where data is cut to fit power decay on
#' lower/upper tail. Set tail.p[1]=0 to exclude lower tail fit; tail.p[2]=1 to exclude upper tail fit.
#' @param transform.info List containing transformation information, returned from \code{ecdfHT}
#' @param x.min Number describing cut-off of lower tail
#' @param x.max Number describing cut-off of upper tail
#' @param weights 'none' to do unweighted regression or 'var' to use weighted regression on tail with weights proportional to variance of quantile
#' @param add.to.plot Boolean indicating whether or not the interpolation is plotted
#' @param ... Optional parameters passed to plot routines, e.g. col='red'
#' @export
#' @return An object of class 'ecdfHT.fit' specifying the interpolation. The fields in the value are:
#' \describe{
#' \item{scale.q}{vector of length 3, copied from the input argument}
#' \item{scale.x}{vector of length 3, the quantiles from the data corresponding to scale.q}
#' \item{xsort}{vector of the sorted, unique data values}
#' \item{ecdf}{nonstandard empirical cdf, see details}
#' \item{xx}{transformed x values: xx[i]=h(xsort[i])}
#' \item{yy}{transformed p values: yy[i]=g(ecdf[i])}
#' \item{cdf.spline}{monotonic spline function used to compute the cdf}
#' \item{inf.cdf.spline}{monotonic spline function used to compute the inverse of the cdf}
#' \item{tail.p}{vector of length 2; probabilities saying where the lower and upper tails begin.
#'     Note these are generally not the exact values of input variable p, rather they are the closest values to those found in ecdf}
#' \item{tail.x}{vector of length 2; x values where the lower and upper tails begin}
#' \item{tail.c}{vector of length 2; tail constants for lower and upper powerlaw fit}
#' \item{tails.slope}{vector of length 2; slope of tails on transformed plot}
#' \item{tail.alpha}{vector of length 2; exponents for lower and upper power law fit}
#' \item{tail.m}{integer vector of length 2; indices in xsort where tails begin }
#' \item{weights}{copy of input variable weights}
#' }
#'
#' @examples
#' x <- rcauchy( 1000 )
#' a <- ecdfHT( x )
#' fit <- ecdfHT.fit( c(.1,.9), a, col='red' )
#' str(fit)
#'
ecdfHT.fit <- function( p, transform.info, x.min=NA, x.max=NA, add.to.plot=TRUE, weights="var", ... ) {

tail.info <- ecdfHT.fit.tails( p, transform.info, weights=weights, add.to.plot=add.to.plot, ... )

n <- length(transform.info$xsort)
# on lower tail back off to guarantee monotonic fit
if ( tail.info$tail.p[1] > 0 ) {
  m1 <- tail.info$tail.m[1]
} else {
  m1 <- 1
}

# on upper tail move right to guarantee monotonic fit
if ( tail.info$tail.p[2] < 1 ) {
  m2 <- tail.info$tail.m[2]
} else {
  m2 <- n
}

# compute the x and y values for the spline fit in the restricted range
xx <- transform.info$xsort[m1:m2]
yy <- transform.info$ecdf[m1:m2]

# create a new object to contain old ...
new <- transform.info

# ... and new information
new$cdf.spline <- approxfun( xx, yy, method="linear",  rule=2  )
new$inv.cdf.spline <- approxfun( yy, xx, method="linear", rule=2 )
new$tail.p <- tail.info$tail.p
new$tail.x  <- tail.info$tail.x
new$tail.c <- tail.info$tail.c
new$tail.alpha <- tail.info$tail.alpha
new$tail.slope <- tail.info$tail.slope
new$tail.m <- c(m1,m2)
new$weights <- weights
class(new) <- "ecdfHT.fit"
return(new ) }
########################################################################
#' @rdname ecdfHT.fit
#' @export
#'
ecdfHT.fit.tails <- function( p, transform.info, weights, add.to.plot=TRUE, ... ) {
#  fit power law model to the lower and upper tails of data.

stopifnot( (weights == "var" | weights == "none"), p[1] < p[2], length(transform.info$xsort) > 3 )
q <- transform.info$scale.q
xx <- transform.info$xx
yy <- transform.info$yy
ecdf <- transform.info$ecdf
xsort <- transform.info$xsort
n <- length(xx)

tail.p <- c(0,1)
tail.x <- c(NA,NA)
tail.slope <- c(NA,NA)
tail.alpha <- c(NA,NA)
tail.c <- c(NA,NA)
m1 <- NA; m2 <- NA
if (p[1] > 0) { # compute a power law on lower/left side
  if (p[1] > q[1]) stop("p[1] > scale.q[1]: lower tail approximation should not be computed on non-logarithmic transformed x values")
  m1 <- findInterval(p[1],ecdf)
  if (m1 < 2) { m1 <- 2 }  # guarantee at least 2 points in regression
  tail.p[1] <- ecdf[m1]
  tail.x[1] <- xsort[m1]
  if( tail.x[1] >= 0 ) stop("error: lower tail cut-off must be negative")

  lower.xx <- xx[1:m1]-xx[m1]
  lower.yy <- yy[1:m1]-yy[m1]
  if (weights=="none") {
    lower.fit <- lm( lower.yy ~ lower.xx - 1 )
  } else {
    w <- (ecdf[1:m1]*(1-ecdf[1:m1]))
    lower.fit <- lm( lower.yy ~ lower.xx - 1, weights=w )
  }
  tail.slope[1] <- lower.fit$coefficients
  tail.alpha[1] <- tail.slope[1]/q[1]
  tail.c[1] <- ecdf[m1]*abs(tail.x[1])^tail.alpha[1]
  if (add.to.plot) {
    y1 <- yy[m1]+tail.slope[1]*(xx[1]-xx[m1])
    lines( c(xx[1],xx[m1]), c(y1,yy[m1]), ... )
    points( xx[m1], yy[m1] , pch="+", cex=3, ... )
  }
}

if (p[2] < 1) { # compute a power law on upper/right side
  if (p[2] < q[3]) stop("p[2] < scale.q[3]: upper tail approximation should not be computed on a non-logarithmic transformed x values")
  m2 <- 1+findInterval(p[2],ecdf,rightmost.closed=TRUE)
  if (m2 >= n) { m2 <- n-1 }  # guarantee at least 2 points in regression
  tail.p[2] <- ecdf[m2]
  tail.x[2] <- xsort[m2]
  if( tail.x[2] <= 0 ) stop("error: upper tail cut-off must be positive")

  upper.xx <- xx[m2:n]-xx[m2]
  upper.yy <- yy[m2:n]-yy[m2]
  if (weights=="none") {
    upper.fit <- lm( upper.yy ~ upper.xx - 1)
  } else {
    w <- (ecdf[m2:n]*(1-ecdf[m2:n]))
    upper.fit <- lm( upper.yy ~ upper.xx - 1, weights=w )
  }
  tail.slope[2] <- upper.fit$coefficients
  tail.alpha[2] <- tail.slope[2]/(1-q[3])
  tail.c[2] <- (1-ecdf[m2])*tail.x[2]^tail.alpha[2]
  if (add.to.plot) {
    y2 <- yy[m2]+tail.slope[2]*(xx[n]-xx[m2])
    lines( c(xx[m2],xx[n]), c(yy[m2],y2), ... )
    points( xx[m2], yy[m2], pch="+", cex=3,... )
  }
}
fit <- list( tail.p=tail.p, tail.x=tail.x, tail.alpha=tail.alpha, tail.slope=tail.slope,
            tail.c=tail.c, tail.m=c(m1,m2) )
invisible( fit ) }
########################################################################
#' Compute cdf, pdf, quantiles, and simulate from a fitted distribution
#'
#' Use the semi-parametric fit calculated by \code{ecdfHT.fit} to evaluate the cdf F(x), pdf f(x), quantiles and simulate.
#' @param x A vector of numbers
#' @param ecdfHT.fit An object returned by \code{ecdfHT.fit} describing the interpolation.
#' @export
#' @details \code{pecdfHT} computes the cdf of the semi-parametric fit to the data.
#' \code{decdfHT} computes the pdf of the semi-parametric fit to the data.  This is likely very irregular and not of much value except on the tails, where the pdf calculation is computed analytically.
#' \code{qecdfHT} computes quantiles.
#' \code{recdfHT} simulates from a semi-parameteric distribution.
#' @return \code{pecdfHT} computes the cdf, \code{decdfHT} computes the pdf, \code{qecdfHT} computes the quantiles (inverse of the cdf), \code{recdfHT} simulates from the distribution.
#' @examples
#' x <- rcauchy(1000)
#' a <- ecdfHT( x, show.plot=FALSE )
#' fit <- ecdfHT.fit( c(.1,.9), a, add.to.plot=FALSE )
#' pecdfHT( -3:3, fit )
#' decdfHT( -3:3, fit )
#' qecdfHT( seq(.1,.9,.1), fit )
#' recdfHT( 10, fit )
#'
pecdfHT <- function( x, ecdfHT.fit ) {
# evaluate the cdf at locations x using spline fit to ecdf in middle
# and power law fit on tails

stopifnot( class(ecdfHT.fit) == "ecdfHT.fit" )

y <- rep( NA, length(x) )
if( is.null(ecdfHT.fit$cdf.spline) ) {
  warning("cannot compute cdf: ecdfHT.fit must be used to setup approximation")
  return(y)
}

# in the center of the data, use spline fit
i <- which( (x >= ecdfHT.fit$tail.x[1]) & (x <= ecdfHT.fit$tail.x[2]) )
y[i] <- ecdfHT.fit$cdf.spline( x[i] )


# lower/left tail power law approximation
i <- which( x < ecdfHT.fit$tail.x[1] )
if (length(i) > 0) {
  y[i] <- ecdfHT.fit$tail.c[1]/abs(x[i])^(ecdfHT.fit$tail.alpha[1])
}

# upper/right tail power law approximation
i <- which( x > ecdfHT.fit$tail.x[2] )
if (length(i) > 0) {
  y[i] <- 1.0 - ecdfHT.fit$tail.c[2] / x[i]^(ecdfHT.fit$tail.alpha[2])
}

return(y) }
########################################################################
#' @rdname pecdfHT
#' @export
decdfHT <- function( x, ecdfHT.fit ) {
# approximate the pdf from the ecdfHT fit
# This is very irregular in the central region where it is piecewise constant, but is
# a smooth power law fit on tails.

stopifnot( class(ecdfHT.fit) == "ecdfHT.fit" )
y <- rep(NA,length(x))
if( is.null(ecdfHT.fit$cdf.spline) ) {
  warning("cannot compute cdf: ecdfHT.fit must be used to setup approximation")
  return(y)
}

# in the middle of the data, compute the derivative of the cdf F by a difference quotient
# since F is piecewise linear, this is an exact calculation
i <- which( (x >= ecdfHT.fit$tail.x[1]) & (x <= ecdfHT.fit$tail.x[2]) )
j <- findInterval( x[i], ecdfHT.fit$xsort, rightmost.closed=TRUE, all.inside=TRUE, left.open=FALSE )
u <- ecdfHT.fit$xsort
v <- ecdfHT.fit$ecdf
y[i] <- (v[j+1]-v[j])/(u[j+1]-u[j])


# lower/left tail
i <- which( x < ecdfHT.fit$tail.x[1] )
if (length(i) > 0) {
  y[i] <- ecdfHT.fit$tail.alpha[1]*ecdfHT.fit$tail.c[1]/abs(x[i])^(ecdfHT.fit$tail.alpha[1]+1)
}

# upper/right tail
i <- which( x > ecdfHT.fit$tail.x[2] )
if (length(i) > 0) {
  y[i] <-  ecdfHT.fit$tail.alpha[2] * ecdfHT.fit$tail.c[2]/x[i]^(ecdfHT.fit$tail.alpha[2]+1)
}

return(y) }
########################################################################
#' @rdname pecdfHT
#' @param p Vector of probabilites
#' @export
#'
qecdfHT <- function( p, ecdfHT.fit ) {
# compute the quantiles from a semi-parametric fit of the data

stopifnot( class(ecdfHT.fit) == "ecdfHT.fit" )

y <- rep( NA, length(p) )
if( is.null(ecdfHT.fit$inv.cdf.spline) ) {
  warning("cannot compute quantiles: ecdfHT.fit must be used to setup approximation")
  return(y)
}

# in the center of the data, use spline fit
i <- which( (p >= ecdfHT.fit$tail.p[1]) & (p <= ecdfHT.fit$tail.p[2]) )
y[i] <- ecdfHT.fit$inv.cdf.spline( p[i] )

# lower/left tail power law approximation
i <- which( p < ecdfHT.fit$tail.p[1] )
if (length(i) > 0) {
  y[i] <- - (ecdfHT.fit$tail.c[1]/p[i])^(1/ecdfHT.fit$tail.alpha[1])
}

# upper/right tail power law approximation
i <- which( p > ecdfHT.fit$tail.p[2] )
if (length(i) > 0) {
  y[i] <- (ecdfHT.fit$tail.c[2]/(1.0-p[i]))^(1/ecdfHT.fit$tail.alpha[2])
}
return(y)}
########################################################################
#' @rdname pecdfHT
#' @param n Number of values to simulate
#' @export
#'
recdfHT <- function( n, ecdfHT.fit ) {
# simulate from the semi-parametric fit of the data

stopifnot( class(ecdfHT.fit) == "ecdfHT.fit" )
u <- runif(n)
x <- qecdfHT( u, ecdfHT.fit )
return(x)}
########################################################################
########################################################################
########################################################################
#' Multivariate extensions of transformed empirical cdf plot
#'
#' Transform multivariate data and plot using the ideas from the univariate plot.
#' @param x Matrix of data of size (n x d)
#' @param q0 quantile of radii transformation
#' @param p.norm Power used in computing L^p norm
#' @param radii.upper.tail.p probability used as cutuoff to tail fit; set to 1 to suppress upper tail fit
#' @param zscale Vector of length 2, value of aspect ratio for the z axis when d=2 and the two 3d plots are drawn
#' @param scale.q matrix of sixe (3 x d), probabilities used to determine the scaling and centering for each component
#' @param show.axes.labels Boolean value, determines if axes are labeled or not
#' @param ... Optional graphical parameters, e.g. col='red'
#' @export
#'
#' @details \code{ecdfHT.multivar} gives a quick graphical look at a d dimensional data set.
#' It produces two plots: the first is a superposition of the univariate \code{ecdfHT} plots
#' for each component; the second plot is an array of plots, showing one plot
#' for each component.
#'
#' \code{ecdfHT.bivar} produces two plots of a bivariate data set.
#' The first one has three subplots:
#' a scatter plot of the data, a transformed scatterplot of the data, and a univariate
#' \code{ecdfHT} plot of the radii of the data.
#' For the second and third subplot, ???? then g(y[,1]) is plotted against g(y[,2]) to get the second plot.
#' For the third plot, compute radius r[i]= l_p norm of shifted and scaled data.  These radii are
#' plotted in a univariate, one-sided \code{ecdfHT} plot.
#'
#' The second plot produced is a 3-dimensonal plot.  It takes the first two subplots just described and adds a
#' third dimension by looking at an ecdf for the radii. Thus the height of a point is low if the point is
#' near the center, and increases as points move away.  The first subplot shows points
#' (x[i,1],x[i,2],ecdf of r[i]).  The second plot transforms all three components:
#' it shows (h1(x[i,1]),hs(x[i,2]),g(ecdf of r[i]), where h1(.) and h2(.) are scaled versions of
#' h(.) from the univariate \code{ecdfHT} plot, and g(.) is as in the univariate plot.
#' See the vignette for more detail.
#'
#'
#' @return \code{ecdfHT.multivar} draws several plots, returns a list (invisibly) with fields:
#' \describe{
#' \item{x}{input (n x d) matrix of data}
#' \item{x.prime}{(n x d) matrix of centered and shifted version of x}
#' \item{y}{(n x d) matrix of transformed x}
#' \item{p.norm}{what p-norm to use; p.norm=2 is Euclidean norm}
#' \item{scale.q}{copy of input argument}
#' \item{radii}{vector of length n, p-norm of the rows of x.prime}
#' \item{q0}{copy of input value}
#' \item{r0}{q0-th quantile of the radii}
#' \item{univariate.ecdfHT}{list of length d, with j-th entry the object returned by \code{ecdfHT} for the j-th column of x}
#' \item{radii.ecdfHT}{list returned from ecdfHT( radii, ... )}
#' \item{radii.tail.fit}{object returned from ecdfHT.fit for the radii )}
#' \item{rgl.id}{rgl id of 3d plot(s); can be used to access, change, print 3d plots}
#' \item{radii.prob}{if d=2, this vector gives the empirical cdf of the radii}
#' \item{radii.prob2}{if d=2, this vector gives the transformed empirical cdf of the radii}
#' }
#'
#' \code{ecdfHT.multivar.transform} computes the transformed vectors y, radii, and
#' \code{lp.norm} computes the lp-norm of the rows of x
#'
#' @examples
#' # independent components
#' set.seed(2)
#' x <- matrix( rcauchy(4000), ncol=4 )
#' ecdfHT.multivar( x )
#'
#' # radially symmetric
#' r <- rcauchy(1000)
#' theta <- runif(1000,min=0,max=2*pi)
#' x <- cbind( r*cos(theta), r*sin(theta) )
#' ecdfHT.multivar( x )

ecdfHT.multivar <- function( x, scale.q=matrix(c(0.25,.5,.75),nrow=3,ncol=ncol(x)), q0=0.5,
                             radii.upper.tail.p=0.90, p.norm=2, show.axes.labels=FALSE,
                             zscale=c(500,1),...){

# compute transform for each coordinate and the radii
multivar.obj <- ecdfHT.multivar.transform( x, scale.q, q0=q0, p.norm=p.norm )
d <- ncol(x)

# multiple plots - one for each coordinate
dev.new( noRStudioGD=TRUE  )
k1 <- ceiling(sqrt(d)); k2 <- ceiling(d/k1)
par(mfrow=c(k1,k2),mar=c(1,1,2,1))
for (j in 1:d) {
  temp <- multivar.obj$univariate.ecdfHT[[j]]
  ecdfHT.draw( temp, temp$xsort, temp$ecdf, new.plot=TRUE,...)
  title(paste("component",j))
}

# single plot with all coordinate overlayed together
dev.new( noRStudioGD=TRUE  )
temp <- multivar.obj$univariate.ecdfHT[[1]]
ecdfHT( temp$xsort, temp$scale.q, type='l', show.axes.labels=FALSE,lwd=3)
for (j in 2:d) {
  temp <- multivar.obj$univariate.ecdfHT[[j]]
  ecdfHT.draw( temp, temp$xsort, temp$ecdf, new.plot=FALSE, col=j,lwd=3)
}
legend("topleft", paste("component",1:d),lty=1,col=1:d, lwd=3, inset=0.05)
title("All components")

# radii ecdfHT
dev.new( noRStudioGD=TRUE  )
temp <- multivar.obj$radii.ecdfHT
ecdfHT( temp$xsort, temp$scale.q, ...)
title("Radii")
if( (radii.upper.tail.p > 0) & (radii.upper.tail.p < 1)) {
   temp2 <- ecdfHT.fit( p=c(0,radii.upper.tail.p), transform.info=temp, col='red', lwd=3 )
   cat("tail fit with radii.upper.tail.p=",radii.upper.tail.p,"\n")
} else {
  temp2 <- NULL
}
multivar.obj$radii.tail.fit <- temp2


# pairs and 2d and 3d cases
rgl.id <- NULL
if (d==2) {
  temp <- ecdfHT.2d( multivar.obj, zscale, ... )
  rgl.id <- temp$rgl.id
  multivar.obj$radii.prob <- temp$radii.prob
  multivar.obj$radii.prob2 <- temp$radii.prob2
} else {
  dev.new( noRStudioGD=TRUE  )
  pairs(x, ...)
  dev.new( noRStudioGD=TRUE  )
  pairs(multivar.obj$y, ...)
  if (d==3){
    rgl.id <- c(-1L,-1L)
    rgl.id[1] <- rgl::open3d()
    rgl::points3d( multivar.obj$x, ... )
    rgl::title3d("Original data X")
    rgl.id[2] <- rgl::open3d( )
    rgl::points3d( multivar.obj$y, ... )
    rgl::title3d("Transformed data Y")
  }
}

multivar.obj$rgl.id <- rgl.id
invisible(multivar.obj)}
########################################################################
#' @rdname ecdfHT.multivar
#' @export
ecdfHT.multivar.transform <- function( x, scale.q, q0, p.norm ) {

# check input
stopifnot( is.matrix(x), ncol(x) > 1, q0 > 0, q0 < 1, p.norm > 0, all(scale.q >= 0),
           all(scale.q<= 1) )
n <- nrow(x)
d <- ncol(x)

# compute all the univariate ecdfHTs (without plotting) and store in a list
univariate.ecdfHT <- vector( mode="list", length= d)
for (j in 1:d) {
  univariate.ecdfHT[[j]] <- ecdfHT( x[,j], scale.q[,j], show.plot=FALSE )
}

# compute radii
x.prime <- matrix(0.0,nrow=nrow(x),ncol=ncol(x))
for (j in 1:d) {
  t <- univariate.ecdfHT[[j]]$scale.x
  x.prime[,j] <- (x[,j]-t[2])/(t[3]-t[2])
}
radii <- lp.norm( x.prime, p.norm=p.norm )
radii.ecdfHT <- ecdfHT(radii, c(0,0,q0), show.plot=FALSE)

# compute transformed y
r0 <- quantile(radii,prob=q0)
h <- ecdfHT.h( radii, c(0,0,r0) )
y <- x.prime

colnames(y) <- paste("y[ ,",1:d,"]")
for (i in 1:n) {
  y[i, ] <- (h[i]/radii[i]) * x.prime[i,]
}

new <- list( x=x, x.prime=x.prime, y=y, p.norm=p.norm, scale.q=scale.q, q0=q0, r0=r0,
             univariate.ecdfHT=univariate.ecdfHT, radii=radii, radii.ecdfHT=radii.ecdfHT )
class(new) <- "ecdfHT.multivar"

invisible(new) }
########################################################################
#' @rdname ecdfHT.multivar
#' @export
#' @param multivar.obj An object of class 'ecdfHT.multivar', see details.
#'
ecdfHT.2d <- function( multivar.obj, zscale=c(500,1), ... ) {

stopifnot( class(multivar.obj) == "ecdfHT.multivar", ncol(multivar.obj$x)==2 )

# x-y scatterplot of raw data
dev.new( noRStudioGD=TRUE  )
par(mfrow=c(1,2),mar=c(2,1,3,1))
plot(multivar.obj$x, xlab="",ylab="", ...)
title("Original data")

# x-y scatterplot of transformed x-y data
plot( multivar.obj$y, ,xlab="",ylab="", ...)
title("Transformed data" )
if (multivar.obj$p.norm==2) {   # show circle where transformation starts
  r0 <- multivar.obj$r0
  theta <- seq(0,2*pi,length=201)
  xx <- r0*cos(theta); yy <- r0*sin(theta)
  lines(xx,yy,lwd=3)
}

# show 3-dim plot of radial cdf
prob <- (rank(multivar.obj$radii) - 0.5)/ length(multivar.obj$radii)
rgl.id <- rgl::open3d()
rgl::points3d( cbind(multivar.obj$x, prob), scale=c(1,1,zscale[1]), ... )
ecdfHT.2d.axes( zscale[1] )   # add x,y,z axes


# 3-d plot of radial cdf in nonlinear x-y scale
# nonlinear scale in probability (z axis)
q0 <-  multivar.obj$q0
i <- which(prob > q0)
prob2 <- prob
prob2[i] <- q0 - log( (1-prob[i])/(1-q0) )
rgl.id[2] <- rgl::open3d()
rgl::points3d( cbind(multivar.obj$y, prob2), ... )
ecdfHT.2d.axes( zscale[2] ) # add x,y,z axes

invisible(list( radii.prob=prob, radii.prob2=prob2, rgl.id=rgl.id ) ) }
########################################################################
#' @rdname ecdfHT.multivar
#' @export
#'
ecdfHT.2d.axes <- function( zscale ){

# add x,y,z axes centered at (0,0,0)
bbox <- rgl::par3d("bbox")
xmax <- max(abs(bbox[1:2]))
rgl::lines3d( c(-xmax,xmax), c(0,0), c(0,0) )
ymax <- max(abs(bbox[3:4]))
rgl::lines3d( c(0,0), c(-ymax,ymax), c(0,0) )
zmax <- max(abs(bbox[5:6]))
rgl::lines3d( c(0,0), c(0,0), c(-.1,1)*zmax )
rgl::aspect3d( 1, 1, zscale*bbox[6]/max(xmax,ymax) )     }
########################################################################
#' @rdname ecdfHT.multivar
#' @export
lp.norm <- function( x, p.norm ){
n <- nrow(x)
r <- rep(0.0,n)
for (i in 1:n) { r[i] <- (sum(x[i,]^p.norm))^(1/p.norm) }
return(r) }
########################################################################
