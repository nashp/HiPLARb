AC_INIT(HiPLARb, version-0.1)

#AC_PROG_CC

MAGMABLAS=''
MAGMA=''
COREBLASSO=''
QUARKSO=''
PLASMASO=''
MAGMAlib=''
BIT64=''
PLASMA=''
PLASMA_HOME=''
MAGMA_HOME=''
CUDA_HOME=''

AC_ARG_WITH([lapack],
    [AC_HELP_STRING([--with-lapack=LIBS],
        [LAPACK shared library @<:@-llapack@:>@])],
   [LIBS="${LIBS} $withval"] [lflag=1], [LIBS="${LIBS} -llapack"] [lflag=0] 
)

if test $lflag -eq 0; then
	  AC_MSG_ERROR([The user should link a version of lapack preferably an optimised version such as ATLAS, OpenBLAS, Intel MKL or AMD ACML. Please see the Install guide for more details.])
fi

#################### PLASMA Libraries ####################
AC_MSG_CHECKING("if PLASMA libraries exist")
AC_ARG_WITH([plasma-lib],
    [AC_HELP_STRING([--with-plasma-lib=DIR],
        [full path to PLASMA shared libraries @<:@@:>@])],
    [PLASMA_HOME="$withval"], AC_MSG_WARN("No PLASMA option given.") 

)
if test ! -z "${PLASMA_HOME}" ; then
	AC_CHECK_FILE([${PLASMA_HOME}/lib/libplasma.so], [PLASMASO="yes"], [PLASMASO=''])
	AC_CHECK_FILE([${PLASMA_HOME}/lib/libcoreblas.so], [COREBLASSO="yes"], [PLASMASO=''])
	AC_CHECK_FILE([${PLASMA_HOME}/lib/libquark.so], [QUARKSO="yes"], [QUARKSO=''])

	if test -z "${PLASMASO}" || test -z "${COREBLASSO}" || test -z "${QUARKSO}" ; then
		  AC_MSG_RESULT("Cannot find dynamic versions of PLASMA.")
	else
		  AC_MSG_RESULT("using ${PLASMA_HOME}/lib for PLASMA_LIBS libs")
		  PLASMA="yes"
		  CFLAGS="${CFLAGS} -I$withval/include -DHIPLAR_WITH_PLASMA"
    	  PLIBS="-lplasma -lcoreblas -lplasma -lcoreblasqw -lquark"
		  LD_FLAGS="${LD_FLAGS} -L${PLASMA_HOME}/lib"
	fi
fi
################# CUDA Lib   ########################

AC_ARG_WITH([cuda-home],
	[AC_HELP_STRING([--with-cuda-home=DIR],
		[full path to CUDA installation @<:@/usr/local/cuda@:>@])],
	[CUDA_HOME="$withval"], AC_MSG_WARN("No CUDA option given.")
)


#################### MAGMA Libraries ####################

AC_ARG_WITH([magma-lib],
	[AC_HELP_STRING([--with-magma-lib=DIR],
		[full path to MAGMA shared libraries @<:@@:>@])],
	[MAGMA_HOME="$withval"], AC_MSG_WARN("No MAGMA option given.")
)

if test ! -z "${MAGMA_HOME}" ; then
	AC_MSG_CHECKING("if dynamic versions of the MAGMA libraries exist")
	AC_CHECK_FILE([${MAGMA_HOME}/lib/libmagma.so],[MAGMA="yes"], [MAGMA=''])

	if test -z "${MAGMA}" ; then
		AC_MSG_RESULT("Cannot find MAGMA libraries")
	else
	  	AC_MSG_RESULT("using ${MAGMA_HOME}/lib for MAGMA_LIBS libs")
	  	MAGMAlib="yes"
	 		LD_FLAGS="${LD_FLAGS} -L${MAGMA_HOME}/lib"
			MLIBS="-lmagma"
		  CFLAGS="${CFLAGS} -I$withval/include -DHIPLAR_WITH_MAGMA -DHAVE_CUBLAS"

		##### Check if CUDA libraries exist ####
  		AC_MSG_CHECKING("if CUDA_HOME exists")
		if test -z "${CUDA_HOME}" ; then
			AC_MSG_ERROR("${CUDA_HOME} does not exist;")
		else
		  	AC_MSG_RESULT("using CUDA_HOME=${CUDA_HOME}")
		fi

		AC_MSG_CHECKING("whether this is the 64 bit linux version of CUDA")
		AC_CHECK_FILE([${CUDA_HOME}/lib64/libcublas.so], [BIT64="yes"])

		if test -z "${BIT64}"; then
			AC_MSG_RESULT("no -- using ${CUDA_HOME}/lib for CUDA libs")
			AC_CHECK_FILE([${CUDA_HOME}/lib/libcublas.so], [LIBCUBLAS="yes"], AC_MSG_ERROR([CUBLAS library not found]))
			CUDA_LIB_DIR="${CUDA_HOME}/lib"
			LD_FLAGS="-L${CUDA_LIB_DIR}"
			CFLAGS="${CFLAGS} -I${CUDA_HOME}/include"
		else
		  	AC_MSG_RESULT("yes -- using ${CUDA_HOME}/lib64 for CUDA libs")
		  	CUDA_LIB_DIR="${CUDA_HOME}/lib64"
		  	LD_FLAGS="${LD_FLAGS} -L${CUDA_LIB_DIR}"
				CFLAGS="${CFLAGS} -I${CUDA_HOME}/include"
		fi
		CLIBS="-lcudart -lcublas -lmagma"
	fi
fi

if test ! -z "${MAGMAlib}" || test ! -z "${PLASMA}" ; then
	  AC_MSG_RESULT("Found one of MAGMA and PLASMA")
else
	 AC_MSG_ERROR([Cannot find MAGMA or PLASMA Libraries. Please install MAGMA, PLASMA or both. Please see the install guide for more details.]) 
fi

LIBS="${LIBS} ${MLIBS} ${CLIBS} ${PLIBS} -lstdc++"

: ${R_HOME=`R RHOME`}
AC_ARG_WITH([r-home],
    [AC_HELP_STRING([--with-r-home=DIR],
        [full path to R installation @<:@@:>@])],
    [R_HOME=$withval]
)
AC_MSG_RESULT([setting R home directory to ${R_HOME}])
AC_CHECK_PROG([check_R], [R], [yes],, [${R_HOME}/bin])
if test -z "${check_R}"; then
    AC_MSG_ERROR([R not found!])
fi

AC_MSG_RESULT([setting R header files to ${R_INCLUDE_DIR}])
CC=`"${R_HOME}/bin/R" CMD config CC`
R_CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
CFLAGS="${CFLAGS} ${R_CPICFLAGS} -I${R_INCLUDE_DIR}"

AC_SUBST([CFLAGS])
AC_SUBST([LD_FLAGS])
AC_SUBST([LIBS])

AC_OUTPUT(src/Makevars)
