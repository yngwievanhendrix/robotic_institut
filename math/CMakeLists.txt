CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(c_libs.math)

IF(NOT INSTALL_DIR)
	FIND_PATH(INSTALL_DIR root.dir 
		${CMAKE_CURRENT_SOURCE_DIR}/.. 
		${CMAKE_CURRENT_SOURCE_DIR}/../.. 
		${CMAKE_CURRENT_SOURCE_DIR}/../../.. 
		${CMAKE_CURRENT_SOURCE_DIR}/../../../.. 
		${CMAKE_CURRENT_SOURCE_DIR}/../../../../.. 
		${CMAKE_CURRENT_SOURCE_DIR}/../../../../../..)
ENDIF(NOT INSTALL_DIR)

SET(LIBRARY_OUTPUT_PATH ${INSTALL_DIR}/lib CACHE PATH "Library output path")

OPTION(BUILD_TEST "Build test program" OFF)
SET(BUILD_TEST ${BUILD_TEST})

SET(ROBMath_TEST
	${CMAKE_CURRENT_SOURCE_DIR}/test/maintest.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/TScal_algorithms.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/TScal_algorithms.h
	${CMAKE_CURRENT_SOURCE_DIR}/test/TSmatrix.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/TSmatrix.h
	${CMAKE_CURRENT_SOURCE_DIR}/test/TShorn.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/test/TShorn.h
)

SET(ROBMath_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/src/wavelets.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/cal_algorithms.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Horn.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/distortion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/matrix.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/vectorops.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/ecodes.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/runningRMS.cpp
)

SET(ROBMath_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/include/wavelets.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/dataobject.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/cal_algorithms.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/Horn.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/matrix.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/vectorops.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/ecodes.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/cal_reply_codes.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/trc_reply_codes.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/rob_reply_codes.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/math.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/distortion.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/runningRMS.h
)

INCLUDE_DIRECTORIES(
	${INSTALL_DIR}/gsl/include
	${INSTALL_DIR}
)

LINK_DIRECTORIES(
  /usr/lib
  /usr/lib64
  /usr/local/lib
	${INSTALL_DIR}/gsl/lib
	${INSTALL_DIR}/lib
)

ADD_LIBRARY(Math
	${ROBMath_HEADERS}
	${ROBMath_SRC}
)

IF(BUILD_TEST MATCHES "^ON$")
	ADD_EXECUTABLE(MathTest
		${ROBMath_TEST}
	)
	TARGET_LINK_LIBRARIES(MathTest
		Math
		cpptest
	)
	LINK_DIRECTORIES(
		${INSTALL_DIR}/cpptest
	)
ENDIF(BUILD_TEST MATCHES "^ON$")


IF(UNIX)
	TARGET_LINK_LIBRARIES(Math
	  gsl
	  gslcblas
	)
	IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" AND NOT BUILD_i386 MATCHES "OFF")
		ADD_LIBRARY(Math-i386
			${ROBMath_Heads}
			${ROBMath_SRC}
		)
		SET_TARGET_PROPERTIES(Math-i386 PROPERTIES COMPILE_FLAGS "-m32 -fPIC")
		TARGET_LINK_LIBRARIES(Math-i386
		  gsl
		  gslcblas
		)
	ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" AND NOT BUILD_i386 MATCHES "OFF")
	SET_TARGET_PROPERTIES( Math PROPERTIES COMPILE_FLAGS "-fPIC" )
ELSE(UNIX)
	IF(CMAKE_GENERATOR MATCHES "Win64")
		TARGET_LINK_LIBRARIES(Math
		  optimized gsl-x64
		  optimized cblas-x64
		  debug gsl-x64_d
		  debug cblas-x64_d
		)
	ELSE(CMAKE_GENERATOR MATCHES "Win64")
		TARGET_LINK_LIBRARIES(Math
		  optimized gsl
		  optimized cblas
		  debug gsl_d
		  debug cblas_d
		)
	ENDIF(CMAKE_GENERATOR MATCHES "Win64")
ENDIF(UNIX)	
