// =========================================================================
// File:       testDeformableContactSurfaceXMLImporter.cpp
// Purpose:    Unit test for VTK mesh loader
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledSelfCollisionBVH.h"
#include "tledSelfCollisionBVHCreator.h"
#include "tledSelfCollisionBVHXMLExporter.h"
#include "tledSelfCollisionBVHXMLImporter.h"
#include "tledGreedySelfCollisionBVHUpdater.h"
#include "tledModel.h"
#include "xmlParser.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>

template <class TSurface>
void _TestXMLWriteRead(TSurface &r_surface) {
  typedef tledSelfCollisionBVHImpl<TSurface, tledAABB<2> > __BVH;

  __BVH input(r_surface), output(r_surface);

  {
    tledSelfCollisionBVHXMLExporter<__BVH> exporter;
    tledSelfCollisionBVHCreator<__BVH> builder;
    tledSelfCollisionBVHXMLImporter<__BVH> importer;

    input.SetMargin(1e-4f);
    input.SetUpdater(*new tledGreedySelfCollisionBVHUpdater<__BVH>());
    input.Init(builder);

    exporter.SetInput(input);
    exporter.Export();
    importer.SetRootNode(exporter.GetRootNode());
    importer.SetOuputObject(output);
    importer.Import();

    exporter.GetRootNode().deleteNodeContent();
  }

  tledUnitTestAssert(input.GetMargin() == output.GetMargin());
 
  tledUnitTestAssert(input.GetLeafBVIndices().size() == output.GetLeafBVIndices().size());
  tledUnitTestAssert(std::equal(input.GetLeafBVIndices().begin(), input.GetLeafBVIndices().end(), output.GetLeafBVIndices().begin()));
 
  tledUnitTestAssert(input.GetNonAdjacentGeometryNodes().size() == output.GetNonAdjacentGeometryNodes().size());
  tledUnitTestAssert(std::equal(input.GetNonAdjacentGeometryNodes().begin(), input.GetNonAdjacentGeometryNodes().end(), output.GetNonAdjacentGeometryNodes().begin()));

  tledUnitTestAssert(input.GetGeometryClusterSubtreeRootIndices().size() == output.GetGeometryClusterSubtreeRootIndices().size());
  tledUnitTestAssert(std::equal(input.GetGeometryClusterSubtreeRootIndices().begin(), input.GetGeometryClusterSubtreeRootIndices().end(), output.GetGeometryClusterSubtreeRootIndices().begin()));

  tledUnitTestAssert(input.GetSelfCollisionCandidates().size() == output.GetSelfCollisionCandidates().size());
  tledUnitTestAssert(std::equal(input.GetSelfCollisionCandidates().begin(), input.GetSelfCollisionCandidates().end(), output.GetSelfCollisionCandidates().begin()));

  tledUnitTestAssert(input.GetNumberOfBVs() == output.GetNumberOfBVs());
  for (int bvInd = 0; bvInd < input.GetNumberOfBVs(); bvInd++) {
    typedef typename __BVH::BoundingVolume __BV;

    const __BV &ibv = input.GetBV(bvInd);
    const __BV &obv = output.GetBV(bvInd);

    tledUnitTestAssert(ibv.ParentIndex == obv.ParentIndex);
    if (ibv.PrimitiveIndex >= 0) {
      tledUnitTestAssert(ibv.PrimitiveIndex == obv.PrimitiveIndex);
    } else {
      tledUnitTestAssert(obv.PrimitiveIndex < 0);
      tledUnitTestAssert(std::equal(ibv.ChildIndices, ibv.ChildIndices + __BV::NumberOfChildBVs, obv.ChildIndices));
    }

    tledUnitTestAssert(ibv.UpdateCounter == obv.UpdateCounter);
    tledUnitTestAssert(std::fabs(ibv.VolinoAngle - obv.VolinoAngle) <= std::fabs(ibv.VolinoAngle)*1e-3);
    if (ibv.VolinoAngle >= 0 && ibv.VolinoAngle <= __BV::GetVolinoThresholdAngle()) {
      tledUnitTestAssert(std::fabs(ibv.SubtreeMinH - obv.SubtreeMinH) <= std::fabs(ibv.SubtreeMinH)*1e-3);
      for (int c = 0; c < 3; c++) {
	tledUnitTestAssert(std::fabs(ibv.VolinoAxis[c] - obv.VolinoAxis[c]) <= std::fabs(ibv.VolinoAxis[c])*1e-3);
      }
    } else {
      tledUnitTestAssert(std::isnan(obv.SubtreeMinH));
      for (int c = 0; c < 3; c++) {
	tledUnitTestAssert(std::isnan(obv.VolinoAxis[c]));
      }
    }

    for (int c = 0; c < 3; c++) for (int b = 0; b < 2; b++) {
      tledUnitTestAssert(std::fabs(ibv.Bounds[c][b] - obv.Bounds[c][b]) <= std::fabs(ibv.Bounds[c][b])*1e-4);
    }
  }
}

static void _TestXMLImport(const std::string &inSimFile) {
  const tledModel model(inSimFile.c_str());
  const tledMesh &mesh = *model.GetMesh();

  if (mesh.GetElType() == std::string("T4") || mesh.GetElType() == std::string("T4ANP")) {
    typedef tledDeformableContactSurfaceT3CPU __Surf;

    __Surf surface(mesh);
    surface.Init();
    _TestXMLWriteRead(surface);
  } else {
    tledFatalError("Incompatible contact surface.");
  }
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestXMLImport(tledUnitTest::GetMeshPath("ball_ball.xml"));

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
