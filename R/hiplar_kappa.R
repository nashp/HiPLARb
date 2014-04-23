#
#  HiPLAR - High Performance Linear Algebra in R
#
#  hiplar_kappa.R is based on src/library/base/R/kappa.R
#  Modifications for the HiPLAR project copyright
#      (C) 2012-2013  Vendel Szeremi, HiPLAR Team
#  HiPLAR is distributed under the same license as R.
#
#  File src/library/base/R/kappa.R
#      Part of the R package, http://www.R-project.org
#      Copyright (C) 1998 B. D. Ripley
#      Copyright (C) 1998-2012 The R Core Team
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


norm <- function(x, type = c("O", "I", "F", "M", "2")) {
    if(identical("2", type)) {
	svd(x, nu=0L, nv=0L)$d[1L]
	## *faster* at least on some platforms {but possibly less accurate}:
	##sqrt(eigen(crossprod(x), symmetric=TRUE, only.values=TRUE)$values[1L])
    } else
	.Call("hiplar_dlange", x, type, PACKAGE="HiPLARb")
}

rcond <- function(x, norm = c("O","I","1"), triangular = FALSE, ...) {
    norm <- match.arg(norm)
    stopifnot(is.matrix(x))
    if({d <- dim(x); d[1L] != d[2L]})## non-square matrix -- use QR
        return(rcond(qr.R(qr(if(d[1L] < d[2L]) t(x) else x)), norm=norm, ...))

    ## x = square matrix :
    if(is.complex(x)) {
        if(triangular)
            .Call("La_ztrcon", x, norm, PACKAGE="base")
        else .Call("hiplar_zgecon", x, norm, PACKAGE="HiPLARb")
    }
    else {
        storage.mode(x) <- "double"
        if(triangular)
            .Call("La_dtrcon", x, norm, PACKAGE="base")
        else .Call("hiplar_dgecon", x, norm, PACKAGE="HiPLARb")
    }
}
