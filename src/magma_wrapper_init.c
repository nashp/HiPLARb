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


#include <stdlib.h>
#include <stdio.h>
#include <R.h>

#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#include <cuda.h>
CUdevice  dev;
CUcontext context;
#endif


int magma_wrapper_init() {
#ifdef HIPLAR_WITH_MAGMA

	magma_init();
/*    if(CUDA_SUCCESS != cuInit(0)) {
        R_ShowMessage("HiPLAR: cuInit failed\n" );
        exit(-1);
    }
*/
/*    if(CUDA_SUCCESS != cuDeviceGet(&dev, 0)) {
        R_ShowMessage("HiPLAR: cuDeviceGet failed\n");
        exit(-1);
    }
*/
 /*   if(CUDA_SUCCESS != cuCtxCreate(&context, 0, dev)) {
        R_ShowMessage("HiPLAR: cuCtxCreate failed");
        exit(-1);
    }*/

    if(CUBLAS_STATUS_SUCCESS != cublasInit()) {
        R_ShowMessage("HiPLAR: cublasInit failed");
        exit(-1);
    }

    magma_print_devices();

#endif
	return 0;
}


int magma_wrapper_deinit() {
#ifdef HIPLAR_WITH_MAGMA

    //cuCtxDetach( context );
    cublasShutdown();

#endif
    return 0;
}
