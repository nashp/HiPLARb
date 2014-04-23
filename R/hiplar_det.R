#
#  HiPLAR - High Performance Linear Algebra in R
#
#  hiplar_det.R is based on src/library/base/R/det.R
#  Modifications for the HiPLAR project copyright
#      (C) 2012-2013  Vendel Szeremi, HiPLAR Team
#  HiPLAR is distributed under the same license as R.
#
#  File src/library/base/R/det.R
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


determinant.matrix <- function(x, logarithm = TRUE, ...)
{
    if ((n <- ncol(x)) != nrow(x))
        stop("'x' must be a square matrix")
    if (n < 1L)
	return(structure(list(modulus =
			      structure(if(logarithm) 0 else 1,
					logarithm = logarithm),
			      sign = 1L),
			 class = "det"))
    if (is.complex(x))
        stop("determinant not currently defined for complex matrices")
    storage.mode(x) <- "double"
    .Call("hiplar_det_ge_real", x, logarithm, PACKAGE = "HiPLARb")
}
