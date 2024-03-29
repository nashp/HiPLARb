#
#  HiPLAR - High Performance Linear Algebra in R
#
#  hiplar_backsolve.R is based on src/library/base/R/backsolve.R
#  Modifications for the HiPLAR project copyright
#      (C) 2012-2013  Vendel Szeremi, HiPLAR Team
#  HiPLAR is distributed under the same license as R.
#
#
#  File src/library/base/R/backsolve.R
#      Part of the R package, http://www.R-project.org
#      Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


backsolve <- function(r, x, k=ncol(r), upper.tri = TRUE, transpose = FALSE)
{
    r <- as.matrix(r) # nr  x  k
    storage.mode(r) <- "double"
    x.mat <- is.matrix(x)
    if(!x.mat) x <- as.matrix(x) # k  x	nb
    storage.mode(x) <- "double"
    k <- as.integer(k)
    ldb <- nrow(x)
    if(is.na(k) || k <= 0 || k > ldb || k > ncol(r))
        stop("invalid argument values in 'backsolve'")
    nb <- as.integer(ncol(x))
    if(is.na(nb)) stop("invalid value of ncol(x)")
    upper.tri <- as.logical(upper.tri)
    transpose <- as.logical(transpose)

#    z <- .Call("hiplar_bakslv2", r, x, ldb, nb, upper.tri, transpose, PACKAGE = "HiPLARb")
#    zx <- if (k < ldb) z[seq_len(k),,drop=FALSE] else z
#    if(x.mat) return(zx) else return(drop(zx))

    job <- as.integer(upper.tri + 10L*transpose)
    z <- .C("hiplar_bakslv",
	    t  = r, ldt= nrow(r), n  = k,
	    b  = x, ldb, nb = nb,
	    x  = matrix(0, ldb, nb),
	    job = job,
	    info = integer(1L),
	    DUP = FALSE, PACKAGE = "HiPLARb")[c("x","info")]

    if(z$info)
	stop(gettextf("singular matrix in 'backsolve'. First zero in diagonal [%d]", z$info), domain = NA)
    zx <- if (k < ldb) z$x[seq_len(k),,drop=FALSE] else z$x
    if(x.mat) zx else drop(zx)

}
