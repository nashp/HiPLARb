/*
 * HiPLAR - High Performance Linear Algebra in R
 *
 * Copyright (C) 2012-2013  Vendel Szeremi, HiPLAR Team
 *
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *          
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *          
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/ 
 */         


#include <string.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include "hiplar_at.h"
#include "hiplar_dbg.h"


int hiplar_library = HIPLAR_USE_AUTO;

int magma_interface = MAGMA_CPU_INTERFACE;

int xover_zgeqrf       = 0;
int xover_dgeqrf       = 0;
int xover_chol         = 0;
int xover_chol2inv     = 0;
int xover_zgesv        = 0;
int xover_dgesv        = 0;
int xover_dlange       = 0;
int xover_bakslv       = 0;
int xover_svd          = 0;
int xover_svd_cmplx    = 0;
int xover_rs           = 0;
int xover_rs_cmplx     = 0;
int xover_rg           = 0;
int xover_rg_cmplx     = 0;
int xover_det_ge_real  = 0;
int xover_dgecon       = 0;
int xover_zgecon       = 0;
int xover_matprod      = 0;

int maxmagma_zgeqrf      = 10000;
int maxmagma_dgeqrf      = 10000;
int maxmagma_chol        = 10000;
int maxmagma_chol2inv    = 10000;
int maxmagma_zgesv       = 10000;
int maxmagma_dgesv       = 10000;
int maxmagma_dlange      = 10000;
int maxmagma_bakslv      = 10000;
int maxmagma_svd         = 10000;
int maxmagma_svd_cmplx   = 10000;
int maxmagma_rs          = 10000;
int maxmagma_rs_cmplx    = 10000;
int maxmagma_rg          = 10000;
int maxmagma_rg_cmplx    = 10000;
int maxmagma_det_ge_real = 10000;
int maxmagma_dgecon      = 10000;
int maxmagma_zgecon      = 10000;
int maxmagma_matprod     = 10000;


SEXP hiplarbSet(SEXP var, SEXP val) {

	const char *pVar = CHAR(STRING_ELT(var,0));
    int tmpVal = asInteger(val);
		
	if(strcmp(pVar, "hiplar_library") == 0) {
        switch (tmpVal) {
            case HIPLAR_USE_PLASMA:
                hiplar_library = HIPLAR_USE_PLASMA;
                break;
            case HIPLAR_USE_MAGMA:
                hiplar_library = HIPLAR_USE_MAGMA;
                break;
            case HIPLAR_USE_AUTO:
                hiplar_library = HIPLAR_USE_AUTO;
                break;
            default:
                hiplar_library = HIPLAR_USE_AUTO;
        }

    } else if (strcmp(pVar, "magma_interface") == 0) {
        magma_interface = tmpVal;
        if ((magma_interface != MAGMA_CPU_INTERFACE) &&
            (magma_interface != MAGMA_GPU_INTERFACE)) {
            magma_interface = MAGMA_CPU_INTERFACE;
        }

	} else if(strcmp(pVar, "xover_zgeqrf") == 0) {
		xover_zgeqrf = tmpVal;
	} else if(strcmp(pVar, "xover_dgeqrf") == 0) {
		xover_dgeqrf = tmpVal;
	} else if(strcmp(pVar, "xover_chol") == 0) {
		xover_chol = tmpVal;
	} else if(strcmp(pVar, "xover_chol2inv") == 0) {
		xover_chol2inv = tmpVal;
	} else if(strcmp(pVar, "xover_zgesv") == 0) {
		xover_zgesv = tmpVal;
	} else if(strcmp(pVar, "xover_dgesv") == 0) {
		xover_dgesv = tmpVal;
	} else if(strcmp(pVar, "xover_dlange") == 0) {
		xover_dlange = tmpVal;
	} else if(strcmp(pVar, "xover_bakslv") == 0) {
		xover_bakslv = tmpVal;
	} else if(strcmp(pVar, "xover_svd") == 0) {
		xover_svd = tmpVal;
	} else if(strcmp(pVar, "xover_svd_cmplx") == 0) {
		xover_svd_cmplx = tmpVal;
	} else if(strcmp(pVar, "xover_rs") == 0) {
		xover_rs = tmpVal;
	} else if(strcmp(pVar, "xover_rs_cmplx") == 0) {
		xover_rs_cmplx = tmpVal;
	} else if(strcmp(pVar, "xover_rg") == 0) {
		xover_rg = tmpVal;
	} else if(strcmp(pVar, "xover_rg_cmplx") == 0) {
		xover_rg_cmplx = tmpVal;
	} else if(strcmp(pVar, "xover_det_ge_real") == 0) {
		xover_det_ge_real = tmpVal;
	} else if(strcmp(pVar, "xover_dgecon") == 0) {
		xover_dgecon = tmpVal;
	} else if(strcmp(pVar, "xover_zgecon") == 0) {
		xover_zgecon = tmpVal;
	} else if(strcmp(pVar, "xover_matprod") == 0) {
		xover_matprod = tmpVal;

	} else if(strcmp(pVar, "maxmagma_zgeqrf") == 0) {
		maxmagma_zgeqrf = tmpVal;
	} else if(strcmp(pVar, "maxmagma_dgeqrf") == 0) {
		maxmagma_dgeqrf = tmpVal;
	} else if(strcmp(pVar, "maxmagma_chol") == 0) {
		maxmagma_chol = tmpVal;
	} else if(strcmp(pVar, "maxmagma_chol2inv") == 0) {
		maxmagma_chol2inv = tmpVal;
	} else if(strcmp(pVar, "maxmagma_zgesv") == 0) {
		maxmagma_zgesv = tmpVal;
	} else if(strcmp(pVar, "maxmagma_dgesv") == 0) {
		maxmagma_dgesv = tmpVal;
	} else if(strcmp(pVar, "maxmagma_dlange") == 0) {
		maxmagma_dlange = tmpVal;
	} else if(strcmp(pVar, "maxmagma_bakslv") == 0) {
		maxmagma_bakslv = tmpVal;
	} else if(strcmp(pVar, "maxmagma_svd") == 0) {
		maxmagma_svd = tmpVal;
	} else if(strcmp(pVar, "maxmagma_svd_cmplx") == 0) {
		maxmagma_svd_cmplx = tmpVal;
	} else if(strcmp(pVar, "maxmagma_rs") == 0) {
		maxmagma_rs = tmpVal;
	} else if(strcmp(pVar, "maxmagma_rs_cmplx") == 0) {
		maxmagma_rs_cmplx = tmpVal;
	} else if(strcmp(pVar, "maxmagma_rg") == 0) {
		maxmagma_rg = tmpVal;
	} else if(strcmp(pVar, "maxmagma_rg_cmplx") == 0) {
		maxmagma_rg_cmplx = tmpVal;
	} else if(strcmp(pVar, "maxmagma_det_ge_real") == 0) {
		maxmagma_det_ge_real = tmpVal;
	} else if(strcmp(pVar, "maxmagma_dgecon") == 0) {
		maxmagma_dgecon = tmpVal;
	} else if(strcmp(pVar, "maxmagma_zgecon") == 0) {
		maxmagma_zgecon = tmpVal;
	} else if(strcmp(pVar, "maxmagma_matprod") == 0) {
		maxmagma_matprod = tmpVal;

    } else {
        R_ShowMessage("ERROR: Unknow variable");
    }

    return R_NilValue;

}


SEXP hiplarbGet(SEXP var) {
	const char *pVar = CHAR(STRING_ELT(var,0));
    SEXP val;
		
    val = PROTECT(allocVector(INTSXP, 1));
	INTEGER(val)[0] = 0;

	if(strcmp(pVar, "hiplar_library") == 0) {
        INTEGER(val)[0] = hiplar_library;
    } else if (strcmp(pVar, "magma_interface") == 0) {
        INTEGER(val)[0] =  magma_interface;

	} else if(strcmp(pVar, "xover_zgeqrf") == 0) {
		INTEGER(val)[0] =xover_zgeqrf;
	} else if(strcmp(pVar, "xover_dgeqrf") == 0) {
		INTEGER(val)[0] =xover_dgeqrf;
	} else if(strcmp(pVar, "xover_chol") == 0) {
		INTEGER(val)[0] =xover_chol;
	} else if(strcmp(pVar, "xover_chol2inv") == 0) {
		INTEGER(val)[0] =xover_chol2inv;
	} else if(strcmp(pVar, "xover_zgesv") == 0) {
		INTEGER(val)[0] =xover_zgesv;
	} else if(strcmp(pVar, "xover_dgesv") == 0) {
		INTEGER(val)[0] =xover_dgesv;
	} else if(strcmp(pVar, "xover_dlange") == 0) {
		INTEGER(val)[0] =xover_dlange;
	} else if(strcmp(pVar, "xover_bakslv") == 0) {
		INTEGER(val)[0] =xover_bakslv;
	} else if(strcmp(pVar, "xover_svd") == 0) {
		INTEGER(val)[0] =xover_svd;
	} else if(strcmp(pVar, "xover_svd_cmplx") == 0) {
		INTEGER(val)[0] =xover_svd_cmplx;
	} else if(strcmp(pVar, "xover_rs") == 0) {
		INTEGER(val)[0] =xover_rs;
	} else if(strcmp(pVar, "xover_rs_cmplx") == 0) {
		INTEGER(val)[0] =xover_rs_cmplx;
	} else if(strcmp(pVar, "xover_rg") == 0) {
		INTEGER(val)[0] =xover_rg;
	} else if(strcmp(pVar, "xover_rg_cmplx") == 0) {
		INTEGER(val)[0] =xover_rg_cmplx;
	} else if(strcmp(pVar, "xover_det_ge_real") == 0) {
		INTEGER(val)[0] =xover_det_ge_real;
	} else if(strcmp(pVar, "xover_dgecon") == 0) {
		INTEGER(val)[0] =xover_dgecon;
	} else if(strcmp(pVar, "xover_zgecon") == 0) {
		INTEGER(val)[0] =xover_zgecon;
	} else if(strcmp(pVar, "xover_matprod") == 0) {
		INTEGER(val)[0] =xover_matprod;

	} else if(strcmp(pVar, "maxmagma_zgeqrf") == 0) {
		INTEGER(val)[0] =maxmagma_zgeqrf;
	} else if(strcmp(pVar, "maxmagma_dgeqrf") == 0) {
		INTEGER(val)[0] =maxmagma_dgeqrf;
	} else if(strcmp(pVar, "maxmagma_chol") == 0) {
		INTEGER(val)[0] =maxmagma_chol;
	} else if(strcmp(pVar, "maxmagma_chol2inv") == 0) {
		INTEGER(val)[0] =maxmagma_chol2inv;
	} else if(strcmp(pVar, "maxmagma_zgesv") == 0) {
		INTEGER(val)[0] =maxmagma_zgesv;
	} else if(strcmp(pVar, "maxmagma_dgesv") == 0) {
		INTEGER(val)[0] =maxmagma_dgesv;
	} else if(strcmp(pVar, "maxmagma_dlange") == 0) {
		INTEGER(val)[0] =maxmagma_dlange;
	} else if(strcmp(pVar, "maxmagma_bakslv") == 0) {
		INTEGER(val)[0] =maxmagma_bakslv;
	} else if(strcmp(pVar, "maxmagma_svd") == 0) {
		INTEGER(val)[0] =maxmagma_svd;
	} else if(strcmp(pVar, "maxmagma_svd_cmplx") == 0) {
		INTEGER(val)[0] =maxmagma_svd_cmplx;
	} else if(strcmp(pVar, "maxmagma_rs") == 0) {
		INTEGER(val)[0] =maxmagma_rs;
	} else if(strcmp(pVar, "maxmagma_rs_cmplx") == 0) {
		INTEGER(val)[0] =maxmagma_rs_cmplx;
	} else if(strcmp(pVar, "maxmagma_rg") == 0) {
		INTEGER(val)[0] =maxmagma_rg;
	} else if(strcmp(pVar, "maxmagma_rg_cmplx") == 0) {
		INTEGER(val)[0] =maxmagma_rg_cmplx;
	} else if(strcmp(pVar, "maxmagma_det_ge_real") == 0) {
		INTEGER(val)[0] =maxmagma_det_ge_real;
	} else if(strcmp(pVar, "maxmagma_dgecon") == 0) {
		INTEGER(val)[0] =maxmagma_dgecon;
	} else if(strcmp(pVar, "maxmagma_zgecon") == 0) {
		INTEGER(val)[0] =maxmagma_zgecon;
	} else if(strcmp(pVar, "maxmagma_matprod") == 0) {
		INTEGER(val)[0] = maxmagma_matprod;

    } else {
        R_ShowMessage("ERROR: Unknow variable");
        return R_NilValue;
    }

	UNPROTECT(1);
	return val;

}


SEXP hiplarShow() {

   Rprintf("HiPLAR library: ");
   switch (hiplar_library) {
        case HIPLAR_USE_PLASMA:
            Rprintf("PLASMA\n");
            break;
        case HIPLAR_USE_MAGMA:
            Rprintf("MAGMA\n");
            break;
        case HIPLAR_USE_AUTO:
            Rprintf("auto\n");
            break;
    }

    Rprintf("xover_zgeqrf                = %5d\n", xover_zgeqrf);
    Rprintf("xover_dgeqrf                = %5d\n", xover_dgeqrf);
    Rprintf("xover_chol                  = %5d\n", xover_chol);
    Rprintf("xover_chol2inv              = %5d\n", xover_chol2inv);
    Rprintf("xover_zgesv                 = %5d\n", xover_zgesv);
    Rprintf("xover_dgesv                 = %5d\n", xover_dgesv);
    Rprintf("xover_dlange                = %5d\n", xover_dlange);
    Rprintf("xover_bakslv                = %5d\n", xover_bakslv);
    Rprintf("xover_svd                   = %5d\n", xover_svd);
    Rprintf("xover_svd_cmplx             = %5d\n", xover_svd_cmplx);
    Rprintf("xover_rs                    = %5d\n", xover_rs);
    Rprintf("xover_rs_cmplx              = %5d\n", xover_rs_cmplx);
    Rprintf("xover_rg                    = %5d\n", xover_rg);
    Rprintf("xover_rg_cmplx              = %5d\n", xover_rg_cmplx);
    Rprintf("xover_det_ge_real           = %5d\n", xover_det_ge_real);
    Rprintf("xover_dgecon                = %5d\n", xover_dgecon);
    Rprintf("xover_zgecon                = %5d\n", xover_zgecon);
    Rprintf("xover_matprod               = %5d\n", xover_matprod);

	Rprintf("maxmagma_zgeqrf             = %5d\n", maxmagma_zgeqrf);
	Rprintf("maxmagma_dgeqrf             = %5d\n", maxmagma_dgeqrf);
	Rprintf("maxmagma_chol               = %5d\n", maxmagma_chol);
	Rprintf("maxmagma_chol2inv           = %5d\n", maxmagma_chol2inv);
	Rprintf("maxmagma_zgesv              = %5d\n", maxmagma_zgesv);
	Rprintf("maxmagma_dgesv              = %5d\n", maxmagma_dgesv);
	Rprintf("maxmagma_dlange             = %5d\n", maxmagma_dlange);
	Rprintf("maxmagma_bakslv             = %5d\n", maxmagma_bakslv);
	Rprintf("maxmagma_svd                = %5d\n", maxmagma_svd);
	Rprintf("maxmagma_svd_cmplx          = %5d\n", maxmagma_svd_cmplx);
	Rprintf("maxmagma_rs                 = %5d\n", maxmagma_rs);
	Rprintf("maxmagma_rs_cmplx           = %5d\n", maxmagma_rs_cmplx);
	Rprintf("maxmagma_rg                 = %5d\n", maxmagma_rg);
	Rprintf("maxmagma_rg_cmplx           = %5d\n", maxmagma_rg_cmplx);
	Rprintf("maxmagma_det_ge_real        = %5d\n", maxmagma_det_ge_real);
	Rprintf("maxmagma_dgecon             = %5d\n", maxmagma_dgecon);
	Rprintf("maxmagma_zgecon             = %5d\n", maxmagma_zgecon);
	Rprintf("maxmagma_matprod            = %5d\n", maxmagma_matprod);

    return R_NilValue;

}


SEXP hiplarWithPLASMA() {

#ifdef HIPLAR_WITH_PLASMA
    return ScalarLogical(TRUE);
#else
    return ScalarLogical(FALSE);
#endif

}


SEXP hiplarWithMAGMA() {

#ifdef HIPLAR_WITH_MAGMA
    return ScalarLogical(TRUE);
#else
    return ScalarLogical(FALSE);
#endif

}


