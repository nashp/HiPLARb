#
#  HiPLAR - High Performance Linear Algebra in R
#
#  hiplar_qr.R is based on src/library/base/R/qr.R
#  Modifications for the HiPLAR project copyright
#      (C) 2012-2013  Vendel Szeremi, HiPLAR Team
#  HiPLAR is distributed under the same license as R.
#
#  File src/library/base/R/qr.R
#      Part of the R package, http://www.R-project.org
#      Copyright (C) 1995-2012 The R Core Team
#
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



#qr.default <- function(x, tol = 1e-07, LAPACK = FALSE, ...)
hiplar_qr <- function(x, tol = 1e-07, LAPACK = FALSE, ...)
{
    x <- as.matrix(x)
    if(is.complex(x)) {
        res <- .Call("hiplar_zgeqrf", x, PACKAGE = "HiPLARb")
	} else {
        if(!is.double(x)) storage.mode(x) <- "double"
        res <- .Call("hiplar_dgeqrf", x, PACKAGE = "HiPLARb")
	}

    if(!is.null(cn <- colnames(x)))
        colnames(res$qr) <- cn[res$pivot]

     class(res) <- "qr"
     return(res)

}


qr.coef <- function(qr, y)
{

    if( !is.qr(qr) )
        stop("first argument must be a QR decomposition")                   

    ap <- attr(qr, "usePLASMA")
    am <- attr(qr, "useMAGMA")
    if ((!is.null(ap) && is.logical(ap) && ap) || (!is.null(am) && is.logical(am) && am)) {

        n <- as.integer(nrow(qr$qr))
        if(is.na(n)) stop("invalid nrow(qr$qr)")
        p <- as.integer(ncol(qr$qr))
        if(is.na(p)) stop("invalid ncol(qr$qr)")

        k <- as.integer(qr$rank)

        im <- is.matrix(y)
        if (!im) y <- as.matrix(y)

        ny <- as.integer(ncol(y))
        if(is.na(ny)) stop("invalid ncol(y)")

        if (p == 0L) return( if (im) matrix(0, p, ny) else numeric() )
            ix <- if ( p > n ) c(seq_len(n), rep(NA, p - n)) else seq_len(p)

        if(is.complex(qr$qr)) {

	        if(!is.complex(y)) y[] <- as.complex(y)
            coef <- matrix(NA_complex_, nrow = p, ncol = ny)
	        coef[qr$pivot,] <-
                .Call("hiplar_coef_cmplx", qr, y, PACKAGE = "HiPLARb")[ix, ]
	        return(if(im) coef else c(coef))

        } else {
            if(!is.double(y)) storage.mode(y) <- "double"
            coef <- matrix(NA_real_, nrow = p, ncol = ny)
            coef[qr$pivot,] <-
                .Call("hiplar_coef_real", qr, y, PACKAGE = "HiPLARb")[ix,]
            return(if(im) coef else c(coef))
		}
	} else {
		return(base::qr.old_coef(qr, y))
	}
}


qr.qy <- function(qr, y)
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")

    ap <- attr(qr, "usePLASMA")
    am <- attr(qr, "useMAGMA")
    if ((!is.null(ap) && is.logical(ap) && ap) || (!is.null(am) && is.logical(am) && am)) {

        if(is.complex(qr$qr)) {
            y <- as.matrix(y)
            if(!is.complex(y)) y[] <- as.complex(y)
            return(.Call("hiplar_qr_qy_cmplx", qr, y, 0, PACKAGE = "HiPLARb"))
        } else {
            return(.Call("hiplar_qr_qy_real", qr, as.matrix(y), 0, PACKAGE = "HiPLARb"))
## storage.mode(y) <- "double"
		}
	} else {
		return(base::qr.old_qy(qr, y))
	}
}


qr.qty <- function(qr, y)
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")

    ap <- attr(qr, "usePLASMA")
    am <- attr(qr, "useMAGMA")
    if ((!is.null(ap) && is.logical(ap) && ap) || (!is.null(am) && is.logical(am) && am)) {

        if(is.complex(qr$qr)){
            y <- as.matrix(y)
            if(!is.complex(y)) y[] <- as.complex(y)
            return(.Call("hiplar_qy_cmplx", qr, y, 1, PACKAGE = "HiPLARb"))
        } else {
            return(.Call("hiplar_qr_qy_real", qr, as.matrix(y), 1, PACKAGE = "HiPLARb"))
		}

	} else {
		return(base::qr.old_qty(qr, y))
	}
}

