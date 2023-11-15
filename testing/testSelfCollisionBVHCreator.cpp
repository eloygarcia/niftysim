// =========================================================================
// File:       testSelfCollisionBVHCreator.cpp
// Purpose:    tledSelfCollisionBVHCreator unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#undef _GPU_

#include "tledUnitTest.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledNarrowConeSelfCollisionBVH.h"
#include "tledSelfCollisionBVHCreator.h"
#include "tledGreedySelfCollisionBVHUpdater.h"
#include "tledModel.h"
#include "tledOBB.h"

template <class TBVH>
static void _GetClusterMaxDepth(int &r_numClusters, int &r_clusterDepth, const TBVH &bvh, const int bvInd, const int myDepth) {
  typedef typename TBVH::BoundingVolume __BV;

  const __BV &bv = bvh.GetBV(bvInd);

  if (!bv.HasDisconnectedGeometry()) {
    r_numClusters += 1;
    r_clusterDepth = std::max(myDepth, r_clusterDepth);
  } else {
    for (int const *pc_c = bv.ChildIndices; pc_c < bv.ChildIndices + __BV::NumberOfChildBVs; pc_c++) if (*pc_c >= 0) {
	_GetClusterMaxDepth(r_numClusters, r_clusterDepth, bvh, *pc_c, myDepth + 1);
      }
  }
}

template <class TBVH>
static void _CheckTopPartBalance(TBVH &r_bvh, const int refDepth, const int numClusters) {
  int numFoundClusters = 0, clusterMaxDepth = 0;

  _GetClusterMaxDepth<TBVH>(numFoundClusters, clusterMaxDepth, r_bvh, 0, 0);
  tledUnitTestAssert(numFoundClusters == numClusters);
  tledUnitTestAssert(clusterMaxDepth == refDepth);
}

template <class TBVH>
static int _GetNumReachableLeafs(const TBVH &bvh, const int bvInd) {
  typedef typename TBVH::BoundingVolume __BV;

  const __BV &bv = bvh.GetBV(bvInd);

  int numLeafs = 0;

  if (bv.PrimitiveIndex >= 0) {
    numLeafs += 1;
  } else {
    for (int const *pc_c = bv.ChildIndices; pc_c < bv.ChildIndices + __BV::NumberOfChildBVs; pc_c++) if (*pc_c >= 0) {
	numLeafs += _GetNumReachableLeafs(bvh, *pc_c);
      }
  }  

  return numLeafs;
}

template <class TBVH, class TBVHBuilder>
static void _CheckBuilder(const std::string &xmlPath, const int refClusterDepth, const int numClusters) {
  typedef typename TBVH::ContactMesh __Surface;

  tledModel model(xmlPath.c_str());
  __Surface surface(*model.GetMesh());

  surface.Init();

  {
    TBVH bvh(surface);
    TBVHBuilder builder;

    bvh.SetMargin(0.05f);
    bvh.SetBVMaxDisplacement(0.02f);
    bvh.SetUpdater(*new tledGreedySelfCollisionBVHUpdater<TBVH>());
    bvh.Init(builder);

    tledUnitTestAssert(_GetNumReachableLeafs(bvh, 0) == surface.GetNumberOfFacets());
    _CheckTopPartBalance(bvh, refClusterDepth, numClusters);
  }
}

static tledUnitTestXMLWriter* _CreateMultiMSHSimulationXML(const int numSubModels) {
  static const std::string boxMsh = tledUnitTest::GetMeshPath("box.msh");
  static const std::string sphereMsh = tledUnitTest::GetMeshPath("sphere.msh");

  tledUnitTestXMLWriter *p_xmlWriter = new tledUnitTestXMLWriter();

  p_xmlWriter->StartXML();
  for (int s = 0; s < numSubModels; s++) {
    const std::string meshFile = s%2 == 1? boxMsh : sphereMsh;

    p_xmlWriter->GetFileStream() << "\t<SubModel>\n";
    
    p_xmlWriter->GetFileStream() << "\t\t<MSHMesh Type=\"T4\">\n"
				 << "\t\t\t" << meshFile << std::endl
				 << "\t\t\t<Translation>\n"
				 << "\t\t\t\t" << s*3 + (s >= 2)*50 << " " << 3*s << " 0\n"
				 << "\t\t\t</Translation>\n"
				 << "\t\t</MSHMesh>\n";
    p_xmlWriter->GetFileStream() << "\t\t<ElementSet Size=\"all\">\n"
				 << "\t\t\t<Material Type=\"NH\">\n"
				 << "\t\t\t\t<ElasticParams NumParams=\"2\">\n"
				 << "\t\t\t\t\t100 1000\n"
				 << "\t\t\t\t</ElasticParams>\n"
				 << "\t\t\t</Material>\n"
				 << "\t\t\t0\n"
				 << "\t\t</ElementSet>\n"
				 << "\t</SubModel>\n";
  }

  p_xmlWriter->GetFileStream() << "\t<SystemParams>\n"
			       << "\t\t<TimeStep>1e-3</TimeStep>\n"
			       << "\t\t<TotalTime>1</TotalTime>\n"
			       << "\t\t<DampingCoeff>0.5</DampingCoeff>\n"
			       << "\t\t<Density>10</Density>\n"
			       << "\t</SystemParams>\n";
  assert(!p_xmlWriter->GetFileStream().fail());
  p_xmlWriter->CloseXML();

  return p_xmlWriter;
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  {
    typedef tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceT3CPU, tledOBB<2> > __BVH;
    typedef tledSelfCollisionBVHCreator<__BVH> __Builder;

    tledUnitTestXMLWriter *p_writer = _CreateMultiMSHSimulationXML(4);

    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 2, 4);
    p_writer->CleanUp();
    delete p_writer;

    p_writer = _CreateMultiMSHSimulationXML(3);
    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 2, 3);
    p_writer->CleanUp();
    delete p_writer;    

    p_writer = _CreateMultiMSHSimulationXML(2);
    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 1, 2);
    p_writer->CleanUp();
    delete p_writer;    

    p_writer = _CreateMultiMSHSimulationXML(1);
    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 0, 1);
    p_writer->CleanUp();
    delete p_writer;    
  }

  {
    typedef tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceT3CPU, tledAABB<4> > __BVH;
    typedef tledSelfCollisionBVHCreator<__BVH> __Builder;

    tledUnitTestXMLWriter *p_writer = _CreateMultiMSHSimulationXML(1);

    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 0, 1);
    p_writer->CleanUp();
    delete p_writer;    

    p_writer = _CreateMultiMSHSimulationXML(2);
    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 1, 2);
    p_writer->CleanUp();
    delete p_writer;    

    //Higher-order BVH top-part construction needs more work!
    p_writer = _CreateMultiMSHSimulationXML(3);
    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 1, 3);
    p_writer->CleanUp();
    delete p_writer;    

    // p_writer = _CreateMultiMSHSimulationXML(4);
    // _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 1, 4);
    // p_writer->CleanUp();
    // delete p_writer;    
  }

  {
    typedef tledNarrowConeSelfCollisionBVH<tledDeformableContactSurfaceT3CPU, tledAABB<2> > __BVH;
    typedef tledSelfCollisionBVHCreator<__BVH> __Builder;

    tledUnitTestXMLWriter *p_writer = _CreateMultiMSHSimulationXML(4);

    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 2, 4);
    p_writer->CleanUp();
    delete p_writer;

    p_writer = _CreateMultiMSHSimulationXML(3);
    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 2, 3);
    p_writer->CleanUp();
    delete p_writer;    

    p_writer = _CreateMultiMSHSimulationXML(2);
    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 1, 2);
    p_writer->CleanUp();
    delete p_writer;    

    p_writer = _CreateMultiMSHSimulationXML(1);
    _CheckBuilder<__BVH, __Builder>(p_writer->GetFilePath(), 0, 1);
    p_writer->CleanUp();
    delete p_writer;    
  }

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
