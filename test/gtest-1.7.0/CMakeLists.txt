
INCLUDE_DIRECTORIES(
	.
	include
)

SET(gtest-1.7.0_SRCS
	src/gtest-all.cc
)

#ADD_DEFINITIONS(-DGTEST_HAS_PTHREAD=1)
ADD_DEFINITIONS(-D_VARIADIC_MAX=10)

FILE(GLOB_RECURSE gtest-1.7.0_HDRS "*.h")

ADD_LIBRARY(gtest ${gtest-1.7.0_SRCS} ${gtest-1.7.0_HDRS})

# Modern MSVC versions no longer provide the legacy std::tr1 tuple API that
# GoogleTest 1.7 probes for.  This test suite does not need TR1 tuples, so
# disable that optional integration for gtest and all of its consumers.
IF(MSVC)
	TARGET_COMPILE_DEFINITIONS(gtest PUBLIC GTEST_HAS_TR1_TUPLE=0)
ENDIF()
