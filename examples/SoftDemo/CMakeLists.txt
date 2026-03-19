INCLUDE_DIRECTORIES(
${BULLET_PHYSICS_SOURCE_DIR}/src
)

LINK_LIBRARIES(
 BulletSoftBody BulletDynamics BulletCollision LinearMath
)

IF (WIN32)
	ADD_EXECUTABLE(App_SoftBodyConsole
		SoftBodyConsole.cpp
		${BULLET_PHYSICS_SOURCE_DIR}/build3/bullet.rc
	)
ELSE()
	ADD_EXECUTABLE(App_SoftBodyConsole
		SoftBodyConsole.cpp
	)
ENDIF()
