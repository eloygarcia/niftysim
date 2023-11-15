// =========================================================================
// File:       testSelfCollisionBVH.cpp
// Purpose:    tledSelfCollisionBVH unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledUnitTest.h"
#include "tledSelfCollisionBVH.h"
#include "tledNarrowConeSelfCollisionBVH.h"
#include "tledSelfCollisionBVHCreator.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledMesh.h"
#include "tledVectorArithmetic.h"
#include "tledGreedySelfCollisionBVHUpdater.h"
#ifdef _Visualisation_
#include "tledVTKMeshLoader.h"
#endif

template <class TBVH>
static void _BuildBVHForSurface(TBVH* (&rp_bvh), typename TBVH::ContactMesh* (&rp_surface), const std::string &meshPath, const std::string &type) {
  typedef typename TBVH::ContactMesh __Mesh;
  
  __Mesh *p_surface;
  TBVH *p_bvh = NULL;

  if (meshPath.substr(meshPath.find_last_of(".") + 1) == "msh") {      
    p_surface = new __Mesh(tledUnitTest::LoadMSHMesh(meshPath, type.c_str()));
  } else {
#ifdef _Visualisation_
    tledMesh solidMesh;
    tledVTKMeshLoader loader;

    loader.SetOutputMesh(solidMesh);
    loader.SetFilename(meshPath);
    loader.SetMeshType(type.c_str());
    loader.Read();

    p_surface = new __Mesh(solidMesh);
#else
    tledFatalFeatureNotEnabledError;
#endif
  }
  p_surface->Init();

  {
    tledSelfCollisionBVHCreator<TBVH> builder;
    
    p_bvh = new TBVH(*p_surface);
    p_bvh->SetMargin(1e-3f);
    p_bvh->SetUpdater(*(new tledGreedySelfCollisionBVHUpdater<TBVH>()));
    p_bvh->Init(builder);
  }

  rp_bvh = p_bvh;
  rp_surface = p_surface;
}

template <class TBVH>
static void _TestConeMonotonyRecursive(const TBVH &bvh, const int bvInd) {  
  tledUnitTestAssert(bvInd >= 0);
  if (bvInd == 0) return;
  else {
    using namespace tledVectorArithmetic;
    typedef typename TBVH::BoundingVolume __BV;

    const __BV &bv = bvh.GetBV(bvInd);
    
    tledUnitTestAssert(bv.ParentIndex >= 0);    
    tledUnitTestAssert(!(bv.HasDisconnectedGeometry() && bvh.IsLeaf(bvInd)));

    {
      const __BV &parent = bvh.GetBV(bv.ParentIndex);

      if (bv.HasDisconnectedGeometry()) {
	tledUnitTestAssert(parent.HasDisconnectedGeometry());
      } else {
	tledUnitTestAssert(bv.VolinoAngle >= 0);

	if (bvh.IsLeaf(bvInd)) {
	  tledUnitTestAssert(bv.VolinoAngle < bv.GetVolinoThresholdAngle());
	}

	if (bv.VolinoAngle >= bv.GetVolinoThresholdAngle()) {
	  tledUnitTestAssert(parent.VolinoAngle >= bv.GetVolinoThresholdAngle());
	} else {
	  tledUnitTestAssert(parent.VolinoAngle >= bv.VolinoAngle);
	  if (parent.VolinoAngle < bv.GetVolinoThresholdAngle()) {
	    tledUnitTestAssert(ComputeAngleNormalised(parent.VolinoAxis, bv.VolinoAxis) <= parent.VolinoAngle);
	  }
	}      
      }
    }
  }
}

template <class TBVH>
static void _TestConeMonotony(const TBVH &bvh) {
  for (std::vector<int>::const_iterator ic_leafInd = bvh.GetLeafBVIndices().begin(); ic_leafInd < bvh.GetLeafBVIndices().end(); ic_leafInd++) {
    _TestConeMonotonyRecursive(bvh, *ic_leafInd);
  }
}

template <class TBVH>
static std::vector<int> _GetSubtreePrimitivesRecursive(const TBVH &bvh, const int root) {
  typedef typename TBVH::BoundingVolume __BV;

  const __BV &bv = bvh.GetBV(root);

  if (bvh.IsLeaf(root)) {
    return std::vector<int>(1, bv.PrimitiveIndex);
  } else {
    std::vector<int> pInds;

    for (int const *pc_c = bv.ChildIndices; pc_c < bv.ChildIndices + __BV::NumberOfChildBVs; pc_c++) if (*pc_c >= 0) {
	std::vector<int> cPInds = _GetSubtreePrimitivesRecursive(bvh, *pc_c);

	pInds.insert(pInds.end(), cPInds.begin(), cPInds.end());
      }

    return pInds;
  }
}

template <class TBVH>
static void _TestSubtreeConeContainment(const TBVH &bvh, const int rootInd) {
  using namespace tledVectorArithmetic;
  typedef typename TBVH::BoundingVolume __BV;
  typedef typename TBVH::ContactMesh __Mesh;

  const __BV &bv = bvh.GetBV(rootInd);  
  
  __Mesh &r_surface = const_cast<__Mesh&>(bvh.GetMesh());

  if (!bvh.IsLeaf(rootInd)) {
    if (!bv.HasDisconnectedGeometry() && bv.VolinoAngle < bv.GetVolinoThresholdAngle()) {
      const std::vector<int> prims = _GetSubtreePrimitivesRecursive(bvh, rootInd);
    
      for (std::vector<int>::const_iterator ic_p = prims.begin(); ic_p < prims.end(); ic_p++) {      
	float n[3], alpha;

	r_surface.ComputeNormalisedFacetNormalCached(n, *ic_p);
	alpha = ComputeAngleNormalised(n, bv.VolinoAxis);
	tledUnitTestAssert(alpha <= bv.VolinoAngle);
      }
    } 
  }

  for (int const *pc_c = bv.ChildIndices; pc_c < bv.ChildIndices + __BV::NumberOfChildBVs; pc_c++) if (*pc_c >= 0) {
      _TestSubtreeConeContainment(bvh, *pc_c);
    }
}

template <class TBVH>
static void _TestAll(const std::string &meshPath, const std::string &meshType) {
  typedef typename TBVH::ContactMesh __Mesh;

  TBVH *p_bvh;
  __Mesh *p_surface;

  _BuildBVHForSurface(p_bvh, p_surface, meshPath, meshType);
  _TestConeMonotony(*p_bvh);
  _TestSubtreeConeContainment(*p_bvh, 0);

  delete p_bvh;
  delete p_surface;  
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestAll<tledSelfCollisionBVHImpl<tledDeformableContactSurfaceT3CPU, tledAABB<2> > >(tledUnitTest::GetMeshPath("sphere.msh"), "T4");
  _TestAll<tledSelfCollisionBVHImpl<tledDeformableContactSurfaceT3CPU, tledAABB<2> > >(tledUnitTest::GetMeshPath("c_extruded_refined.msh"), "T4");

  _TestAll<tledSelfCollisionBVHImpl<tledDeformableContactSurfaceT3CPU, tledAABB<4> > >(tledUnitTest::GetMeshPath("sphere.msh"), "T4");
  _TestAll<tledSelfCollisionBVHImpl<tledDeformableContactSurfaceT3CPU, tledAABB<4> > >(tledUnitTest::GetMeshPath("c_extruded_refined.msh"), "T4");

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
