/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2012 Erwin Coumans  http://bulletphysics.org

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


#ifndef BT_BULLET_XML_WORLD_IMPORTER_H
#define BT_BULLET_XML_WORLD_IMPORTER_H

#include "LinearMath/btScalar.h"

class btDynamicsWorld;

namespace tinyxml2
{
class XMLNode;
};

struct btConvexInternalShapeData;
struct btCollisionShapeData;
#ifdef BT_USE_DOUBLE_PRECISION
struct btRigidBodyDoubleData;
struct btTypedConstraintDoubleData;
#define btRigidBodyData btRigidBodyDoubleData
#define btTypedConstraintData2 btTypedConstraintDoubleData
#else
struct btRigidBodyFloatData;
struct btTypedConstraintFloatData;
#define btTypedConstraintData2 btTypedConstraintFloatData
#define btRigidBodyData btRigidBodyFloatData
#endif  //BT_USE_DOUBLE_PRECISION

struct btCompoundShapeChildData;

#include "LinearMath/btAlignedObjectArray.h"
#include "btWorldImporter.h"

class btBulletXmlWorldImporter : public btWorldImporter
{
protected:
	btAlignedObjectArray<btCollisionShapeData*> m_collisionShapeData;
	btAlignedObjectArray<btAlignedObjectArray<btCompoundShapeChildData>*> m_compoundShapeChildDataArrays;
	btAlignedObjectArray<btRigidBodyData*> m_rigidBodyData;
	btAlignedObjectArray<btTypedConstraintData2*> m_constraintData;
	btHashMap<btHashPtr, void*> m_pointerLookup;
	int m_fileVersion;
	bool m_fileOk;

	void auto_serialize_root_level_children(tinyxml2::XMLNode* pParent);
	void auto_serialize(tinyxml2::XMLNode* pParent);

	void deSerializeVector3FloatData(tinyxml2::XMLNode* pParent, btAlignedObjectArray<btVector3FloatData>& vectors);

	void fixupCollisionDataPointers(btCollisionShapeData* shapeData);
	void fixupConstraintData(btTypedConstraintData2* tcd);

	//collision shapes data
	void deSerializeCollisionShapeData(tinyxml2::XMLNode* pParent, btCollisionShapeData* colShapeData);
	void deSerializeConvexInternalShapeData(tinyxml2::XMLNode* pParent);
	void deSerializeStaticPlaneShapeData(tinyxml2::XMLNode* pParent);
	void deSerializeCompoundShapeData(tinyxml2::XMLNode* pParent);
	void deSerializeCompoundShapeChildData(tinyxml2::XMLNode* pParent);
	void deSerializeConvexHullShapeData(tinyxml2::XMLNode* pParent);
	void deSerializeDynamicsWorldData(tinyxml2::XMLNode* parent);

	///bodies
	void deSerializeRigidBodyFloatData(tinyxml2::XMLNode* pParent);

	///constraints
	void deSerializeGeneric6DofConstraintData(tinyxml2::XMLNode* pParent);

public:
	btBulletXmlWorldImporter(btDynamicsWorld* world);

	virtual ~btBulletXmlWorldImporter();

	bool loadFile(const char* fileName);
};

#endif  //BT_BULLET_XML_WORLD_IMPORTER_H
