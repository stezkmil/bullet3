#ifndef BT_DEFORMABLE_VBD_SETTINGS_H
#define BT_DEFORMABLE_VBD_SETTINGS_H
#include "LinearMath/btScalar.h"
struct btDeformableVbdSettings
{
	int iterations = 10;
	int workers = 1;
	btScalar collisionUnitsPerMeter = 1000;
	btScalar radius = .0001, gap = .001, ke = 1e6, kd = 10, frictionEpsilon = .01, volumeFloor = .05;
	btScalar recoveryDistance = .0001;
};
#endif
