/*
 * HiPLAR - High Performance Linear Algebra in R
 *
 * util.c contains code from src/modules/lapack/Lapack.c
 * HiPLAR is distributed under the same license as R.
 *
 * src/modules/lapack/Lapack.c
 *     R : A Computer Language for Statistical Data Analysis
 *     Copyright (C) 2001--2012  The R Core Team.
 *     Copyright (C) 2003--2010  The R Foundation
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
#include <ctype.h> /* for toupper */

#include <Rdefines.h>

/* Lapack condition number approximation: currently only supports _1 or _Inf norm : */
char La_rcond_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error("argument type[1]='%s' must be a character string of string length 1",
	      typstr);
    typup = toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* alias */
    else if (typup != 'O' && typup != 'I')
	error("argument type[1]='%s' must be one of '1','O', or 'I'",
	      typstr);
    return typup; /* 'O' or 'I' */
}
