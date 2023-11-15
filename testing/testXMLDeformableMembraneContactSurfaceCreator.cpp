// =========================================================================
// File:       testXMLDeformableMembraneContactSurfaceCreator.cpp
// Purpose:    tledXMLDeformableMembraneContactSurfaceCreator unit test
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
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledDeformableMembraneContactSurfaceXMLExporter.h"
#include "tledXMLDeformableMembraneContactSurfaceCreator.h"

template <class TSurface>
void _TestWriteRead(TSurface &r_surface, const float t[]) {
  using namespace tledVectorArithmetic;

  TSurface output;

  {
    tledDeformableMembraneContactSurfaceXMLExporter<TSurface> exporter;
    tledXMLDeformableMembraneContactSurfaceCreator<TSurface> importer;

    exporter.SetInput(r_surface);
    exporter.Export();
    importer.SetXMLRoot(exporter.GetRootNode());
    importer.SetOutputMesh(output);
    importer.Create();
    output.SetMembraneFacetThickenesses(t);

    exporter.GetRootNode().deleteNodeContent();
  }  

  tledUnitTestAssert(output.GetNumberOfMembraneNodes() == r_surface.GetNumberOfMembraneNodes());
  tledUnitTestAssert(output.GetNumberOfMembraneElements() == r_surface.GetNumberOfMembraneElements());
  tledUnitTestAssert(output.GetNumberOfMembraneEdges() == r_surface.GetNumberOfMembraneEdges());

  tledUnitTestAssert(output.GetMembraneNodeBaseIndex() == r_surface.GetMembraneNodeBaseIndex());
  tledUnitTestAssert(output.GetMembraneEdgeBaseIndex() == r_surface.GetMembraneEdgeBaseIndex());
  tledUnitTestAssert(output.GetMembraneFacetBaseIndex() == r_surface.GetMembraneFacetBaseIndex());
  
  for (int n = 0; n < r_surface.GetNumberOfNodes(); n++) {
    float d[3];

    tledUnitTestAssert(Norm(Sub(d, r_surface.GetNodeNormalCached(n), output.GetNodeNormalCached(n))) <= 1e-3);
    tledUnitTestAssert(Norm(Sub(d, r_surface.GetNodeCoordinates(n), output.GetNodeCoordinates(n))) <= 1e-3);
    tledUnitTestAssert(Norm(Sub(d, r_surface.GetNodeCoordinates0(n), output.GetNodeCoordinates0(n))) <= 1e-3);
  }

  for (int f = 0; f < r_surface.GetNumberOfFacets(); f++) {
    float d[3];

    tledUnitTestAssert(Norm(Sub(d, r_surface.GetNormalisedOldFacetNormalCached(f), output.GetNormalisedOldFacetNormalCached(f))) <= 1e-3);
  }
}

static void _TestXMLExportImport(const std::string &inSimFile) {
  const tledModel model(inSimFile.c_str());
  const tledMesh &mesh = *model.GetMesh();

  assert(model.GetNumberOfShellElementSets() > 0 && !model.DoShellUseMeshSurface());
  if (mesh.GetElType() == std::string("T4") || mesh.GetElType() == std::string("T4ANP")) {
    typedef tledDeformableMembraneContactSurfaceT3CPU __Surf;
      
    tledSurface *p_membrane = model.GetGenericShellMesh();
    std::vector<float> t(p_membrane->GetNumberOfFacets(), 0.3f);
    __Surf surface(mesh, *p_membrane);
      
    surface.SetMembraneFacetThickenesses(&t.front());
    surface.Init();
    _TestWriteRead(surface, &t.front());

    delete p_membrane;
  } else {
    tledFatalError("Incompatible contact surface.");
  }
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestXMLExportImport(tledUnitTest::GetMeshPath("two_cubes_sheet.xml"));
  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
