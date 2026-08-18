
# These suites use utilities that are provided by Extras (notably
# BulletInverseDynamicsUtils).  Do not create broken test targets when a
# library-only build, such as the Conan recipe, disables BUILD_EXTRAS.
IF(BUILD_BULLET3 AND BUILD_EXTRAS)
	#SUBDIRS( TestBullet3OpenCL )
	SUBDIRS(  InverseDynamics SharedMemory )
ENDIF(BUILD_BULLET3 AND BUILD_EXTRAS)

SUBDIRS(  gtest-1.7.0 collision BulletDynamics )

