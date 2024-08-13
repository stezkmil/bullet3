/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2009 Erwin Coumans  http://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_TRIANGLE_CALLBACK_H
#define BT_TRIANGLE_CALLBACK_H

#include "LinearMath/btVector3.h"

///The btTriangleCallback provides a callback for each overlapping triangle when calling processAllTriangles.
///This callback is called by processAllTriangles for all btConcaveShape derived class, such as  btBvhTriangleMeshShape, btStaticPlaneShape and btHeightfieldTerrainShape.
class btTriangleCallback
{
	bool m_useProcessTriangleEx = false;

public:
	btTriangleCallback()
	{

	}
	btTriangleCallback(bool useProcessTriangleEx) : m_useProcessTriangleEx(useProcessTriangleEx)
	{

	}
	virtual ~btTriangleCallback();
	bool usesProcessTriangleEx() const
	{
		return m_useProcessTriangleEx;
	}
	virtual void processTriangle(btVector3* triangle, int partId, int triangleIndex) = 0;
	virtual void processTriangleEx(btVector3* triangle, int* triangleVertexIndex, int partId, int triangleIndex)
	{
		processTriangle(triangle, partId, triangleIndex);
	}
};

class btInternalTriangleIndexCallback
{
public:
	virtual ~btInternalTriangleIndexCallback();
	virtual void internalProcessTriangleIndex(btVector3* triangle, int partId, int triangleIndex) = 0;
};

#endif  //BT_TRIANGLE_CALLBACK_H
