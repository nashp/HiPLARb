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


hiplarb_mode_plasma <- function() {
    .Call("hiplarbSet", "hiplar_library", 1, PACKAGE="HiPLARb")
}


hiplarb_mode_magma <- function() {
    .Call("hiplarbSet", "hiplar_library", 2, PACKAGE="HiPLARb")
}


hiplarb_mode_auto <- function() {
    .Call("hiplarbSet", "hiplar_library", 3, PACKAGE="HiPLARb")
}


hiplarb_magma_cpu_int <- function() {
    .Call("hiplarbSet", "magma_interface", 0, PACKAGE="HiPLARb")
}


hiplarb_magma_gpu_int <- function() {
    .Call("hiplarbSet", "magma_interface", 1, PACKAGE="HiPLARb")
}


hiplarb_set_param <- function(param) {

    .Call("hiplarbSet", "xover_dgeqrf", param["dgeqrf", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_dgeqrf", param["dgeqrf", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_zgeqrf", param["zgeqrf", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_zgeqrf", param["zgeqrf", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_chol", param["chol", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_chol", param["chol", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_chol2inv", param["chol2inv", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_chol2inv", param["chol2inv", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_zgesv", param["zgesv", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_zgesv", param["zgesv", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_dgesv", param["dgesv", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_dgesv", param["dgesv", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_dlange", param["dlange", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_dlange", param["dlange", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_bakslv", param["bakslv", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_bakslv", param["bakslv", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_svd", param["svd", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_svd", param["svd", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_svd_cmplx", param["svd_cmplx", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_svd_cmplx", param["svd_cmplx", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_rs", param["rs", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_rs", param["rs", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_rs_cmplx", param["rs_cmplx", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_rs_cmplx", param["rs_cmplx", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_rg", param["rg", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_rg", param["rg", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_rg_cmplx", param["rg_cmplx", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_rg_cmplx", param["rg_cmplx", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_det_ge_real", param["det_ge_real", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_det_ge_real", param["det_ge_real", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_dgecon", param["dgecon", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_dgecon", param["dgecon", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_zgecon", param["zgecon", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_zgecon", param["zgecon", "maxmagma"], PACKAGE="HiPLARb")

    .Call("hiplarbSet", "xover_matprod", param["matprod", "xover"], PACKAGE="HiPLARb")
    .Call("hiplarbSet", "maxmagma_matprod", param["matprod", "maxmagma"], PACKAGE="HiPLARb")
}


hiplarb_get_param <- function() {

    hiplarb_param <- data.frame(xover=array(NA,18), maxmagma=array(NA,18))
    row.names(hiplarb_param) <- c(
    	"dgeqrf", "zgeqrf",    "chol",        "chol2inv", 
    	"zgesv",  "dgesv",     "dlange",      "bakslv", 
    	"svd",    "svd_cmplx", "rs",          "rs_cmplx", 
    	"rg",     "rg_cmplx",  "det_ge_real", "dgecon", 
    	"zgecon", "matprod")

    hiplarb_param["dgeqrf", "xover"] <- .Call("hiplarbGet", "xover_dgeqrf", PACKAGE="HiPLARb")
    hiplarb_param["zgeqrf", "xover"] <- .Call("hiplarbGet", "xover_zgeqrf", PACKAGE="HiPLARb")
    hiplarb_param["chol", "xover"] <- .Call("hiplarbGet", "xover_chol", PACKAGE="HiPLARb")
    hiplarb_param["chol2inv", "xover"] <- .Call("hiplarbGet", "xover_chol2inv", PACKAGE="HiPLARb")
    hiplarb_param["zgesv", "xover"] <- .Call("hiplarbGet", "xover_zgesv", PACKAGE="HiPLARb")
    hiplarb_param["dgesv", "xover"] <- .Call("hiplarbGet", "xover_dgesv", PACKAGE="HiPLARb")
    hiplarb_param["dlange", "xover"] <- .Call("hiplarbGet", "xover_dlange", PACKAGE="HiPLARb")
    hiplarb_param["bakslv", "xover"] <- .Call("hiplarbGet", "xover_bakslv", PACKAGE="HiPLARb")
    hiplarb_param["svd", "xover"] <- .Call("hiplarbGet", "xover_svd", PACKAGE="HiPLARb")
    hiplarb_param["svd_cmplx", "xover"] <- .Call("hiplarbGet", "xover_svd_cmplx", PACKAGE="HiPLARb")
    hiplarb_param["rs", "xover"] <- .Call("hiplarbGet", "xover_rs", PACKAGE="HiPLARb")
    hiplarb_param["rs_cmplx", "xover"] <- .Call("hiplarbGet", "xover_rs_cmplx", PACKAGE="HiPLARb")
    hiplarb_param["rg", "xover"] <- .Call("hiplarbGet", "xover_rg", PACKAGE="HiPLARb")
    hiplarb_param["rg_cmplx", "xover"] <- .Call("hiplarbGet", "xover_rg_cmplx", PACKAGE="HiPLARb")
    hiplarb_param["det_ge_real", "xover"] <- .Call("hiplarbGet", "xover_det_ge_real", PACKAGE="HiPLARb")
    hiplarb_param["dgecon", "xover"] <- .Call("hiplarbGet", "xover_dgecon", PACKAGE="HiPLARb")
    hiplarb_param["zgecon", "xover"] <- .Call("hiplarbGet", "xover_zgecon", PACKAGE="HiPLARb")
    hiplarb_param["matprod", "xover"] <- .Call("hiplarbGet", "xover_matprod", PACKAGE="HiPLARb")

    hiplarb_param["dgeqrf", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_dgeqrf", PACKAGE="HiPLARb")
    hiplarb_param["zgeqrf", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_zgeqrf", PACKAGE="HiPLARb")
    hiplarb_param["chol", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_chol", PACKAGE="HiPLARb")
    hiplarb_param["chol2inv", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_chol2inv", PACKAGE="HiPLARb")
    hiplarb_param["zgesv", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_zgesv", PACKAGE="HiPLARb")
    hiplarb_param["dgesv", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_dgesv", PACKAGE="HiPLARb")
    hiplarb_param["dlange", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_dlange", PACKAGE="HiPLARb")
    hiplarb_param["bakslv", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_bakslv", PACKAGE="HiPLARb")
    hiplarb_param["svd", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_svd", PACKAGE="HiPLARb")
    hiplarb_param["svd_cmplx", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_svd_cmplx", PACKAGE="HiPLARb")
    hiplarb_param["rs", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_rs", PACKAGE="HiPLARb")
    hiplarb_param["rs_cmplx", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_rs_cmplx", PACKAGE="HiPLARb")
    hiplarb_param["rg", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_rg", PACKAGE="HiPLARb")
    hiplarb_param["rg_cmplx", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_rg_cmplx", PACKAGE="HiPLARb")
    hiplarb_param["det_ge_real", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_det_ge_real", PACKAGE="HiPLARb")
    hiplarb_param["dgecon", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_dgecon", PACKAGE="HiPLARb")
    hiplarb_param["zgecon", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_zgecon", PACKAGE="HiPLARb")
    hiplarb_param["matprod", "maxmagma"] <- .Call("hiplarbGet", "maxmagma_matprod", PACKAGE="HiPLARb")

    return(hiplarb_param)

}


hiplarb_do_xover_timings <- function(min=10, max=3000, rep=20) {

    if ((.Call("hiplarWithPLASMA", PACKAGE="HiPLARb") == FALSE) ||
        (.Call("hiplarWithMAGMA", PACKAGE="HiPLARb") == FALSE)) {
        print("Only available when compiled with both PLASMA and MAGMA.")
        return()
    }

    print(sprintf("   Started  %s", Sys.time()))

    param <- hiplarb_get_param()
    hiplarb_magma_cpu_int()

    # use PLASMA
    param["dlange", "xover"] <- param["dlange", "maxmagma"]
    param["rs", "xover"] <- param["rs", "maxmagma"]
    param["rs_cmplx", "xover"] <- param["rs_cmplx", "maxmagma"]

    # use MAGMA
    param["svd", "xover"] <- 0
    param["svd_cmplx", "xover"] <- 0
    param["rg", "xover"] <- 0
    param["rg_cmplx", "xover"] <- 0


    x_initial <- min
    if (min < 100) {
        b <- 90
        if (max < 100) { b <- max }
        x_initial <- c(x_initial, seq(min+10, b, 10))
    }
    if ((min < 1000) && (max > 100)) {
        b <- 900
        if (max < 1000) { b <- max }
        if (min < 100) {
            x_initial <- c(x_initial, seq(100, b, 100))
        } else {
            x_initial <- seq(min, b, 100)
        }
    }
    if (max > 1000) {
        if (min < 1000) {
            x_initial <- c(x_initial, seq(1000, max, 500))
        } else {
            x_initial <- seq(min, max, 500)
        }
    }
    if (x_initial[length(x_initial)] != max) {
        x_initial <- c(x_initial, max)
    }


    ## dgesv
    x <- x_initial 

    mt <- rep(0.0, length(x))
    for (j in length(x):1) {
        i <- x[j]
        A <- matrix(data=2.0,nrow=i,ncol=i)
        .Call("plasma_wrapper_dplrnt", A, 51, PACKAGE = "HiPLARb")
        B <- matrix(data=2.0,nrow=i,ncol=1)
        .Call("plasma_wrapper_dplrnt", A, 5673, PACKAGE = "HiPLARb")
        b <- c(B)
            
        hiplarb_mode_plasma()
        t <- replicate(rep, system.time(solve(A,b), gcFirst = TRUE)[3])
        pmean <- mean(t)
        hiplarb_mode_magma()
        t <- replicate(rep, system.time(solve(A,b), gcFirst = TRUE)[3])
        mmean <- mean(t)

        #print(sprintf("size=%4d p=%f m=%f", i, pmean, mmean))
        mt[j] = pmean - mmean

        rm(A)
        rm(B)
        rm(b)
        if ((mt[j] < 0) && (j < length(x))) {
            break
        }
    }

    if (mt[length(x)] < 0) {
        print(sprintf("dgesv    search terminated - try larger size for max argument; using size=%4d ", x[length(x)]))
        i <- x[length(x)]
    } else if (mt[j] > 0) {
        print(sprintf("dgesv    search terminated - try smaller size for min argument; using size=%4d ", x[j]))
        i <- x[j]
    } else {
        if (j == length(x)) { j <- j - 1 }
        inc <- as.integer((x[j+1] - x[j]) / 10)
        if (inc < 1) { inc <- 1 }
        x <- seq(x[j], x[j+1], inc)
        mt <- rep(0.0, length(x))

        for (j in 1:length(x)) {
            i <- x[j]
            A <- matrix(data=2.0,nrow=i,ncol=i)
            .Call("plasma_wrapper_dplrnt", A, 51, PACKAGE = "HiPLARb")
            B <- matrix(data=2.0,nrow=i,ncol=1)
            .Call("plasma_wrapper_dplrnt", A, 5673, PACKAGE = "HiPLARb")
            b <- c(B)
            
            hiplarb_mode_plasma()
            t <- replicate(rep, system.time(solve(A,b), gcFirst = TRUE)[3])
            pmean <- mean(t)
            hiplarb_mode_magma()
            t <- replicate(rep, system.time(solve(A,b), gcFirst = TRUE)[3])
            mmean <- mean(t)

            mt[j] = pmean - mmean
            rm(A)
            rm(B)
            rm(b)
        }
        for (j in length(x):1) {
            if (mt[j] < 0) {
                break
            }
        }
        i <- x[j]
        print(sprintf("dgesv    size=%4d", i))
    }

    param["dgesv", "xover"] <- i
    param["dgecon", "xover"] <- i
    param["det_ge_real", "xover"] <- i


    ## zgesv
    x <- x_initial 

    mt <- rep(0.0, length(x))
    for (j in length(x):1) {
        i <- x[j]

        A <- matrix(data=as.complex(2.0),nrow=i,ncol=i)
        .Call("plasma_wrapper_zplrnt", A, 51, PACKAGE = "HiPLARb")
        B <- matrix(data=as.complex(2.0),nrow=i,ncol=1)
        .Call("plasma_wrapper_zplrnt", A, 5673, PACKAGE = "HiPLARb")
        b <- c(B)

        hiplarb_mode_plasma()
        t <- replicate(rep, system.time(solve(A,b), gcFirst = TRUE)[3])
        pmean <- mean(t)
        hiplarb_mode_magma()
        t <- replicate(rep, system.time(solve(A,b), gcFirst = TRUE)[3])
        mmean <- mean(t)

        #print(sprintf("size=%4d p=%f m=%f", i, pmean, mmean))
        mt[j] = pmean - mmean

        rm(A)
        rm(B)
        rm(b)
        if ((mt[j] < 0) && (j < length(x))) {
            break
        }
    }

    if (mt[length(x)] < 0) {
        print(sprintf("zgesv    search terminated - try larger size for max argument; using size=%4d ", x[length(x)]))
        i <- x[length(x)]
    } else if (mt[j] > 0) {
        print(sprintf("zgesv    search terminated - try smaller size for min argument; using size=%4d ", x[j]))
        i <- x[j]
    } else {
        if (j == length(x)) { j <- j - 1 }
        inc <- as.integer((x[j+1] - x[j]) / 10)
        if (inc < 1) { inc <- 1 }
        x <- seq(x[j], x[j+1], inc)
        mt <- rep(0.0, length(x))

        for (j in 1:length(x)) {
            i <- x[j]

            A <- matrix(data=as.complex(2.0),nrow=i,ncol=i)
            .Call("plasma_wrapper_zplrnt", A, 51, PACKAGE = "HiPLARb")
            B <- matrix(data=as.complex(2.0),nrow=i,ncol=1)
            .Call("plasma_wrapper_zplrnt", A, 5673, PACKAGE = "HiPLARb")
            b <- c(B)

            hiplarb_mode_plasma()
            t <- replicate(rep, system.time(solve(A,b), gcFirst = TRUE)[3])
            pmean <- mean(t)
            hiplarb_mode_magma()
            t <- replicate(rep, system.time(solve(A,b), gcFirst = TRUE)[3])
            mmean <- mean(t)

            mt[j] = pmean - mmean
            rm(A)
            rm(B)
            rm(b)
        }
        for (j in length(x):1) {
            if (mt[j] < 0) {
                break
            }
        }
        i <- x[j]
        print(sprintf("zgesv    size=%4d", i))
    }

    param["zgesv", "xover"] <- i
    param["zgecon", "xover"] <- i


    ## chol
    x <- x_initial 

    mt <- rep(0.0, length(x))
    for (j in length(x):1) {
        i <- x[j]

        B <- matrix(data=2.0,nrow=i,ncol=i)
        .Call("plasma_wrapper_dplgsy", B, 51, PACKAGE = "HiPLARb")

        hiplarb_mode_plasma()
        t <- replicate(rep, system.time(chol(B), gcFirst = TRUE)[3])
        pmean <- mean(t)
        hiplarb_mode_magma()
        t <- replicate(rep, system.time(chol(B), gcFirst = TRUE)[3])
        mmean <- mean(t)

        #print(sprintf("size=%4d p=%f m=%f", i, pmean, mmean))
        mt[j] = pmean - mmean

        rm(B)
        if ((mt[j] < 0) && (j < length(x))) {
            break
        }
    }

    if (mt[length(x)] < 0) {
        print(sprintf("chol     search terminated - try larger size for max argument; using size=%4d ", x[length(x)]))
        i <- x[length(x)]
    } else if (mt[j] > 0) {
        print(sprintf("chol     search terminated - try smaller size for min argument; using size=%4d ", x[j]))
        i <- x[j]
    } else {
        if (j == length(x)) { j <- j - 1 }
        inc <- as.integer((x[j+1] - x[j]) / 10)
        if (inc < 1) { inc <- 1 }
        x <- seq(x[j], x[j+1], inc)
        mt <- rep(0.0, length(x))

        for (j in 1:length(x)) {
            i <- x[j]

            B <- matrix(data=2.0,nrow=i,ncol=i)
            .Call("plasma_wrapper_dplgsy", B, 51, PACKAGE = "HiPLARb")

            hiplarb_mode_plasma()
            t <- replicate(rep, system.time(chol(B), gcFirst = TRUE)[3])
            pmean <- mean(t)
            hiplarb_mode_magma()
            t <- replicate(rep, system.time(chol(B), gcFirst = TRUE)[3])
            mmean <- mean(t)

            mt[j] = pmean - mmean
            rm(B)
        }
        for (j in length(x):1) {
            if (mt[j] < 0) {
                break
            }
        }
        i <- x[j]
        print(sprintf("chol     size=%4d", i))
    } 

    param["chol", "xover"] <- i


    ## chol2inv
    x <- x_initial 

    mt <- rep(0.0, length(x))
    for (j in length(x):1) {
        i <- x[j]

        B <- matrix(data=2.0,nrow=i,ncol=i)
        .Call("plasma_wrapper_dplgsy", B, 51, PACKAGE = "HiPLARb")
        A <- chol(B)

        hiplarb_mode_plasma()
        t <- replicate(rep, system.time(chol2inv(A), gcFirst = TRUE)[3])
        pmean <- mean(t)
        hiplarb_mode_magma()
        t <- replicate(rep, system.time(chol2inv(A), gcFirst = TRUE)[3])
        mmean <- mean(t)

        #print(sprintf("size=%4d p=%f m=%f", i, pmean, mmean))
        mt[j] = pmean - mmean

        rm(A)
        rm(B)
        if ((mt[j] < 0) && (j < length(x))) {
            break
        }
    }

    if (mt[length(x)] < 0) {
        print(sprintf("chol2inv search terminated - try larger size for max argument; using size=%4d ", x[length(x)]))
        i <- x[length(x)]
    } else if (mt[j] > 0) {
        print(sprintf("chol2inv search terminated - try smaller size for min argument; using size=%4d ", x[j]))
        i <- x[j]
    } else {
        if (j == length(x)) { j <- j - 1 }
        inc <- as.integer((x[j+1] - x[j]) / 10)
        if (inc < 1) { inc <- 1 }
        x <- seq(x[j], x[j+1], inc)
        mt <- rep(0.0, length(x))

        for (j in 1:length(x)) {
            i <- x[j]

            B <- matrix(data=2.0,nrow=i,ncol=i)
            .Call("plasma_wrapper_dplgsy", B, 51, PACKAGE = "HiPLARb")
            A <- chol(B)

            hiplarb_mode_plasma()
            t <- replicate(rep, system.time(chol2inv(A), gcFirst = TRUE)[3])
            pmean <- mean(t)
            hiplarb_mode_magma()
            t <- replicate(rep, system.time(chol2inv(A), gcFirst = TRUE)[3])
            mmean <- mean(t)

            mt[j] = pmean - mmean
            rm(A)
            rm(B)
        }
        for (j in length(x):1) {
            if (mt[j] < 0) {
                break
            }
        }
        i <- x[j]
        print(sprintf("chol2inv size=%4d", i))
    }

    param["chol2inv", "xover"] <- i

if(0){
    ## qr: dgeqrf
    x <- x_initial 

    mt <- rep(0.0, length(x))
    for (j in length(x):1) {
        i <- x[j]

        B <- matrix(data=2.0,nrow=i,ncol=i)
        .Call("plasma_wrapper_dplrnt", B, 3456, PACKAGE = "HiPLARb")

        hiplarb_mode_plasma()
#        t <- replicate(rep, system.time(qr(B), gcFirst = TRUE)[3])
        t <- replicate(rep, system.time(hiplar_qr(B), gcFirst = TRUE)[3])
        pmean <- mean(t)
        hiplarb_mode_magma()
#        t <- replicate(rep, system.time(qr(B), gcFirst = TRUE)[3])
        t <- replicate(rep, system.time(hiplar_qr(B), gcFirst = TRUE)[3])
        mmean <- mean(t)

        #print(sprintf("size=%4d p=%f m=%f", i, pmean, mmean))
        mt[j] = pmean - mmean

        rm(B)
        if ((mt[j] < 0) && (j < length(x))) {
            break
        }
    }

    if (mt[length(x)] < 0) {
        print(sprintf("dgeqrf   search terminated - try larger size for max argument; using size=%4d ", x[length(x)]))
        i <- x[length(x)]
    } else if (mt[j] > 0) {
        print(sprintf("dgeqrf   search terminated - try smaller size for min argument; using size=%4d ", x[j]))
        i <- x[j]
    } else {
        if (j == length(x)) { j <- j - 1 }
        inc <- as.integer((x[j+1] - x[j]) / 10)
        if (inc < 1) { inc <- 1 }
        x <- seq(x[j], x[j+1], inc)
        mt <- rep(0.0, length(x))

        for (j in 1:length(x)) {
            i <- x[j]
            B <- matrix(data=2.0,nrow=i,ncol=i)
            .Call("plasma_wrapper_dplrnt", B, 3456, PACKAGE = "HiPLARb")

            hiplarb_mode_plasma()
#        t <- replicate(rep, system.time(qr(B), gcFirst = TRUE)[3])
            t <- replicate(rep, system.time(hiplar_qr(B), gcFirst = TRUE)[3])
            pmean <- mean(t)
            hiplarb_mode_magma()
#        t <- replicate(rep, system.time(qr(B), gcFirst = TRUE)[3])
            t <- replicate(rep, system.time(hiplar_qr(B), gcFirst = TRUE)[3])
            mmean <- mean(t)

            mt[j] = pmean - mmean
            rm(B)
        }
        for (j in length(x):1) {
            if (mt[j] < 0) {
                break
            }
        }
        i <- x[j]
        print(sprintf("dgeqrf   size=%4d", i))
    }

    param["dgeqrf", "xover"] <- i
} else {
    param["dgeqrf", "xover"] <- param["dgeqrf", "maxmagma"]
}

    ## qr: zgeqrf
if (0) {
    x <- x_initial 

    mt <- rep(0.0, length(x))
    for (j in length(x):1) {
        i <- x[j]

        B <- matrix(data=as.complex(2.0),nrow=i,ncol=i)
        .Call("plasma_wrapper_zplrnt", B, 3456, PACKAGE = "HiPLARb")

        hiplarb_mode_plasma()
#        t <- replicate(rep, system.time(qr(B), gcFirst = TRUE)[3])
        t <- replicate(rep, system.time(hiplar_qr(B), gcFirst = TRUE)[3])
        pmean <- mean(t)
        hiplarb_mode_magma()
#        t <- replicate(rep, system.time(qr(B), gcFirst = TRUE)[3])
        t <- replicate(rep, system.time(hiplar_qr(B), gcFirst = TRUE)[3])
        mmean <- mean(t)

        #print(sprintf("size=%4d p=%f m=%f", i, pmean, mmean))
        mt[j] = pmean - mmean

        rm(B)
        if ((mt[j] < 0) && (j < length(x))) {
            break
        }
    }

    if (mt[length(x)] < 0) {
        print(sprintf("bakslv   search terminated - try larger size for max argument; using size=%4d ", x[length(x)]))
        i <- x[length(x)]
    } else if (mt[j] > 0) {
        print(sprintf("bakslv   search terminated - try smaller size for min argument; using size=%4d ", x[j]))
        i <- x[j]
    } else {
        if (j == length(x)) { j <- j - 1 }
        inc <- as.integer((x[j+1] - x[j]) / 10)
        if (inc < 1) { inc <- 1 }
        x <- seq(x[j], x[j+1], inc)
        mt <- rep(0.0, length(x))

        for (j in 1:length(x)) {
            i <- x[j]

            B <- matrix(data=as.complex(2.0),nrow=i,ncol=i)
            .Call("plasma_wrapper_zplrnt", B, 3456, PACKAGE = "HiPLARb")

            hiplarb_mode_plasma()
#        t <- replicate(rep, system.time(qr(B), gcFirst = TRUE)[3])
            t <- replicate(rep, system.time(hiplar_qr(B), gcFirst = TRUE)[3])
            pmean <- mean(t)
            hiplarb_mode_magma()
#        t <- replicate(rep, system.time(qr(B), gcFirst = TRUE)[3])
            t <- replicate(rep, system.time(hiplar_qr(B), gcFirst = TRUE)[3])
            mmean <- mean(t)

            mt[j] = pmean - mmean
            rm(B)
        }
        for (j in length(x):1) {
            if (mt[j] < 0) {
                break
            }
        }
        i <- x[j]
        print(sprintf("zgeqrf   size=%4d", i))
    }

    param["zgeqrf", "xover"] <- i
} else {
    param["zgeqrf", "xover"] <- param["zgeqrf", "maxmagma"]
}

    ## bakslv
    x <- x_initial 

    mt <- rep(0.0, length(x))
    for (j in length(x):1) {
        i <- x[j]

        A <- matrix(data=2.0,nrow=i,ncol=i)
        .Call("plasma_wrapper_dplrnt", A, 453, PACKAGE = "HiPLARb")
        B <- matrix(data=2.0,nrow=i,ncol=1)
        .Call("plasma_wrapper_dplrnt", A, 5673, PACKAGE = "HiPLARb")
        b <- c(B)

        backsolve(A,b)

        hiplarb_mode_plasma()
        t <- replicate(rep, system.time(backsolve(A,b), gcFirst = TRUE)[3])
        pmean <- mean(t)
        hiplarb_mode_magma()
        t <- replicate(rep, system.time(backsolve(A,b), gcFirst = TRUE)[3])
        mmean <- mean(t)

        #print(sprintf("size=%4d p=%f m=%f", i, pmean, mmean))
        mt[j] = pmean - mmean

        rm(A)
        rm(B)
        rm(b)

        if ((mt[j] < 0) && (j < length(x))) {
            break
        }
    }

    if (mt[length(x)] < 0) {
        print(sprintf("bakslv   search terminated - try larger size for max argument; using size=%4d ", x[length(x)]))
        i <- x[length(x)]
    } else if (mt[j] > 0) {
        print(sprintf("bakslv   search terminated - try smaller size for min argument; using size=%4d ", x[j]))
        i <- x[j]
    } else {
        if (j == length(x)) { j <- j - 1 }
        inc <- as.integer((x[j+1] - x[j]) / 10)
        if (inc < 1) { inc <- 1 }
        x <- seq(x[j], x[j+1], inc)
        mt <- rep(0.0, length(x))

        for (j in 1:length(x)) {
            i <- x[j]

            A <- matrix(data=2.0,nrow=i,ncol=i)
            .Call("plasma_wrapper_dplrnt", A, 453, PACKAGE = "HiPLARb")
            B <- matrix(data=2.0,nrow=i,ncol=1)
            .Call("plasma_wrapper_dplrnt", A, 5673, PACKAGE = "HiPLARb")
            b <- c(B)

            hiplarb_mode_plasma()
            t <- replicate(rep, system.time(backsolve(A,b), gcFirst = TRUE)[3])
            pmean <- mean(t)
            hiplarb_mode_magma()
            t <- replicate(rep, system.time(backsolve(A,b), gcFirst = TRUE)[3])
            mmean <- mean(t)

            mt[j] = pmean - mmean

            rm(A)
            rm(B)
            rm(b)
        }

        for (j in length(x):1) {
            if (mt[j] < 0) {
                break
            }
        }
        i <- x[j]
        print(sprintf("bakslv   size=%4d", i))
    }

    param["bakslv", "xover"] <- i


#    ## rs
#    i <- min
#    nf <- TRUE
#    while (nf) {    
#        B <- matrix(data=2.0,nrow=i,ncol=i)
#        .Call("plasma_wrapper_dplgsy", B, 51, PACKAGE = "HiPLARb")
#
#        hiplarb_mode_plasma()
#        t <- replicate(rep, system.time(eigen(B, symmetric=TRUE), gcFirst = TRUE)[3])
#        pmean <- mean(t)
#        hiplarb_mode_magma()
#        t <- replicate(rep, system.time(eigen(B, symmetric=TRUE), gcFirst = TRUE)[3])
#        mmean <- mean(t)
#
#        print(sprintf("rs       size=%4d Plasma=%f Magma=%f", i, pmean, mmean))
#
#        if (pmean > mmean) {
#            nf <- FALSE
#        }
#        rm(B)
#        i <- i + inc
#    }
#    param["rs", "xover"] <- i


    ## matprod
    x <- x_initial 

    mt <- rep(0.0, length(x))
    for (j in length(x):1) {
        i <- x[j]

        A <- matrix(data=2.0,nrow=i,ncol=i)
        .Call("plasma_wrapper_dplrnt", A, 453, PACKAGE = "HiPLARb")
        B <- matrix(data=2.0,nrow=i,ncol=i)
        .Call("plasma_wrapper_dplrnt", B, 5673, PACKAGE = "HiPLARb")

        hiplarb_mode_plasma()
        t <- replicate(rep, system.time(A %*% B, gcFirst = TRUE)[3])
        pmean <- mean(t)
        hiplarb_mode_magma()
        t <- replicate(rep, system.time(A %*% B, gcFirst = TRUE)[3])
        mmean <- mean(t)

        #print(sprintf("size=%4d p=%f m=%f", i, pmean, mmean))
        mt[j] = pmean - mmean

        rm(A)
        rm(B)
        if ((mt[j] < 0) && (j < length(x))) {
            break
        }
    }

    if (mt[length(x)] < 0) {
        print(sprintf("matprod  search terminated - try larger size for max argument; using size=%4d ", x[length(x)]))
        i <- x[length(x)]
    } else if (mt[j] > 0) {
        if (x[j] < 20) {
            print(sprintf("matprod  search terminated; using size=%4d ", x[j]))
        } else {
            print(sprintf("matprod  search terminated - try smaller size for min argument; using size=%4d ", x[j]))
        }
        i <- x[j]
    } else {
        if (j == length(x)) { j <- j - 1 }
        inc <- as.integer((x[j+1] - x[j]) / 10)
        if (inc < 1) { inc <- 1 }
        x <- seq(x[j], x[j+1], inc)
        mt <- rep(0.0, length(x))

        for (j in 1:length(x)) {
            i <- x[j]

            A <- matrix(data=2.0,nrow=i,ncol=i)
            .Call("plasma_wrapper_dplrnt", A, 453, PACKAGE = "HiPLARb")
            B <- matrix(data=2.0,nrow=i,ncol=i)
            .Call("plasma_wrapper_dplrnt", B, 5673, PACKAGE = "HiPLARb")

            hiplarb_mode_plasma()
            t <- replicate(rep, system.time(A %*% B, gcFirst = TRUE)[3])
            pmean <- mean(t)
            hiplarb_mode_magma()
            t <- replicate(rep, system.time(A %*% B, gcFirst = TRUE)[3])
            mmean <- mean(t)

            mt[j] = pmean - mmean
            rm(A)
            rm(B)
        }
        for (j in length(x):1) {
            if (mt[j] < 0) {
                break
            }
        }
        i <- x[j]
        print(sprintf("matprod  size=%4d", i))
    }

    param["matprod", "xover"] <- i

    print(sprintf("   Finished %s", Sys.time()))

    write.table(param, file="~/.hiplarb.csv")
    return(param)
}


