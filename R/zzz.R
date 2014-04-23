#
#  HiPLAR - High Performance Linear Algebra in R
#
#      (C) 2012-2013  Vendel Szeremi, HiPLAR Team
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


.onLoad <- function(libname, pkgname) {
	.C("hiplar_init", PACKAGE="HiPLARb")


################################################################################
# Replace the do_matprod function in the R evaluator.
# This requires src/main/objects.c to be patched:
# attribute_hidden needs to be removed from R_possible_dispatch.
# 
# SEXP attribute_hidden
# R_possible_dispatch(SEXP call, SEXP op, SEXP args, SEXP rho,
#            Rboolean promisedArgs)

    .Call("hiplar_replace")

################################################################################


	if (getRversion() < "2.14.0") {
		aIN <- assignInNamespace
	} else {
		aIN <- function(x, value, ns, pos = -1, envir = as.environment(pos)) {
		    if(missing(ns)) {
		        nm <- attr(envir, "name", exact = TRUE)
		        if(is.null(nm) || substring(nm, 1L, 8L) != "package:")
		            stop("environment specified is not a package")
		        ns <- asNamespace(substring(nm, 9L))
		    } else ns <- asNamespace(ns)

		    if(bindingIsLocked(x, ns)) {
		        unlockBinding(x, ns)
		        assign(x, value, envir = ns, inherits = FALSE)
		        w <- options("warn")
		        on.exit(options(w))
		        options(warn = -1)
		        lockBinding(x, ns)
		    } else {
		        assign(x, value, envir = ns, inherits = FALSE)
		    }
		    if(!isBaseNamespace(ns)) {
		        ## now look for possible copy as a registered S3 method
		        S3 <- getNamespaceInfo(ns, "S3methods")
		        if(!length(S3)) return(invisible(NULL))
		        S3names <- S3[, 3L]
		        if(x %in% S3names) {
		            i <- match(x, S3names)
		            genfun <- get(S3[i, 1L], mode = "function", envir = parent.frame())
		            if(.isMethodsDispatchOn() && methods::is(genfun, "genericFunction"))
		                genfun <- methods::slot(genfun, "default")@methods$ANY
		            defenv <- if (typeof(genfun) == "closure") environment(genfun)
		            else .BaseNamespaceEnv
		            S3Table <- get(".__S3MethodsTable__.", envir = defenv)
		            remappedName <- paste(S3[i, 1L], S3[i, 2L], sep = ".")
		            if(exists(remappedName, envir = S3Table, inherits = FALSE))
		                assign(remappedName, value, S3Table)
		        }
		    }
		    invisible(NULL)
		}
	}

	aIN("qr.old_qy", base::qr.qy, ns="base")
	aIN("qr.qy", qr.qy, ns="base")
	aIN("qr.old_qty", base::qr.qty, ns="base")
	aIN("qr.qty", qr.qty, ns="base")
	aIN("qr.old_coef", base::qr.coef, ns="base")
	aIN("qr.coef", qr.coef, ns="base")

#	aIN("qr.old_default", base::qr.default, ns="base")
#	aIN("qr.default", qr.default, ns="base")

	aIN("chol.old_default", base::chol.default, ns="base")
	aIN("chol.default", chol.default, ns="base")
	aIN("old_chol2inv", base::chol2inv, ns="base")
	aIN("chol2inv", chol2inv, ns="base")

	aIN("solve.old_default", base::solve.default, ns="base")
	aIN("solve.default", solve.default, ns="base")

	aIN("old_backsolve", base::backsolve, ns="base")
	aIN("backsolve", backsolve, ns="base")

	aIN("La.old_svd", base::La.svd, ns="base")
	aIN("La.svd", La.svd, ns="base")

	aIN("old_eigen", base::eigen, ns="base")
	aIN("eigen", eigen, ns="base")

	if (getRversion() >= "2.11.0") {
		aIN("old_norm", base::norm, ns="base")
		aIN("norm", norm, ns="base")
	}

	aIN("determinant.old_matrix", base::determinant.matrix, ns="base")
	aIN("determinant.matrix", determinant.matrix, ns="base")

    aIN("old_rcond", base::rcond, ns="base")
    aIN("rcond", rcond, ns="base")

}


.onAttach <- function(libname, pkgname) {

    if (file.exists("~/.hiplarb.csv")) {
        a <- read.table("~/.hiplarb.csv")
        hiplarb_set_param(a)
    }

}


.onUnload <- function(libname) {
	.C("hiplar_deinit", PACKAGE="HiPLARb")

	assignInNamespace("qr.qy", base::qr.old_qy, ns="base")
	assignInNamespace("qr.qty", base::qr.old_qty, ns="base")
	assignInNamespace("qr.coef", base::qr.old_coef, ns="base")
#	assignInNamespace("qr.default", base::qr.old_default, ns="base")

	assignInNamespace("chol.default", base::chol.old_default, ns="base")
	assignInNamespace("chol2inv", base::old_chol2inv, ns="base")

	assignInNamespace("solve.default", base::solve.old_default, ns="base")

	assignInNamespace("backsolve", base::old_backsolve, ns="base")

	assignInNamespace("La.svd", base::La.old_svd, ns="base")

	assignInNamespace("eigen", base::old_eigen, ns="base")

	if (getRversion() >= "2.11.0") {
		assignInNamespace("norm", base::old_norm, ns="base")
	}

	assignInNamespace("determinant.matrix", base::determinant.old_matrix, ns="base")

    assignInNamespace("rcond", base::old_rcond, ns="base")
}

