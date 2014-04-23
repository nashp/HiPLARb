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


#include "plasma_wrapper_init.h"
#include "magma_wrapper_init.h"

void hiplar_init() {

	plasma_wrapper_init();
    magma_wrapper_init();

}

void hiplar_deinit() {

	plasma_wrapper_deinit();
    magma_wrapper_deinit();

}
